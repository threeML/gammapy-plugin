import numpy as np
from threeML.io.logging import setup_logger
from gammapy.modeling.models import SkyModel, SpectralModel, SpatialModel, TemporalModel
from gammapy.modeling.parameter import Parameter, Parameters
from astromodels.core.model import Model
from astromodels.sources import PointSource, ExtendedSource, Source
from astromodels.functions.function import Function, Function1D, Function2D, Function3D


log = setup_logger(__name__)


class AstromodelConverter:
    """
    Class for analyizing an astromodel model and converting all the individual
    sources such that it can be used with gammapy

    Every Source in the Model will get its own Gammapy skymodel.
    The evaluation happens via the astromodel definition.

    !!! Currently only single source Models supported !!!
    """

    def __init__(self, model: Model, frame: str = None) -> None:
        """
        :param model: the astromodel model
        :type model: astromodel.core.model.Model
        :param frame: geometry frame for gammapy, defaults to ICRS
        :type frame: str
        :return:

        :example:

        """
        assert isinstance(model, Model), "Needs an astromodels Model"
        self._astromodel_model = model
        self._frame = frame

        self._converted_sources = {}
        self._gammapy_models = []
        self._convert_extendend_sources()
        self._convert_point_sources()
        self._create_gammapy_models()

    def _convert_extendend_sources(self) -> None:
        """
        Converts an extended source into a gammapy skymodel
        """
        for (
            source_name,
            source_instance,
        ) in self._astromodel_model.extended_sources.items():
            raise NotImplementedError("currently fails")

            self._converted_sources[source_name] = SourceConverter(
                source_instance, frame=self._frame, converter=self
            )

    def _convert_point_sources(self) -> None:
        """
        Converts a point source into a skymodel
        """
        for (
            source_name,
            source_instance,
        ) in self._astromodel_model.point_sources.items():
            self._converted_sources[source_name] = SourceConverter(
                source_instance, converter=self
            )

    def _create_gammapy_models(self) -> None:
        for name, source in self._converted_sources.items():
            self._gammapy_models.append(source.skymodel)

    def _update_parameters(self) -> None:
        """
        Update all the
        """
        for name, source in self._converted_sources.items():
            for pn, p in self._astromodel_model.free_parameters.items():
                if name in pn:
                    source._update_parameter(pn, p.value)

    @property
    def gammapy_models(self) -> list[SkyModel]:
        """
        Returns all the gammapy skymodels for that model
        """
        return self._gammapy_models


class SourceConverter:
    def __init__(
        self, source: Source, converter: AstromodelConverter, **kwargs
    ) -> None:
        self._frame = kwargs.get("frame", None)
        self._converter = converter
        self._spatial_model = None
        self._spectral_model = None
        self._temporal_model = None
        self._source = source
        self._all_para_names = list(self._source.parameters.keys())

        if isinstance(self._source, PointSource):
            log.debug("Source is a PointSource")
            self._convert_spectral_model()

        elif isinstance(self._source, ExtendedSource):
            log.debug("Source is an ExtendedSource")
            self._convert_spatial_model()
            self._convert_spectral_model()
        else:
            log.error("This astromodels source type is not yet supported.")
            raise NotImplementedError

        self._create_skymodel()

    def _create_parameter_dict(self) -> None:
        """
        Initially creates the parameter dict for that source
        Only the value will be updated during the sampling
        which should be faster
        """
        # TODO thats likely inefficient
        self._parameter_dict = {}
        # This initalizes the parameter dict for an indiviudal source
        for name in self._skymodel.parameters.names:
            self._parameter_dict[name] = {}
            astromodel_para = self._converter._astromodel_model.parameters[name]
            self._parameter_dict[name]["value"] = astromodel_para.value
            self._parameter_dict[name]["unit"] = astromodel_para.unit
            val = np.nan
            if astromodel_para.min_value is not None:
                val = astromodel_para.min_value
            self._parameter_dict[name]["min"] = val
            val = np.nan
            if astromodel_para.max_value is not None:
                val = astromodel_para.max_value
            self._parameter_dict[name]["max"] = val
            self._parameter_dict[name]["frozen"] = not astromodel_para.free
            self._parameter_dict[name]["prior"] = ""

    def _convert_spectral_model(self) -> None:
        """
        Convert the spectral model of the source
        """
        spectral_models = []
        for comp_name, comp in self._source.components.items():
            para_names = []
            for p in self._all_para_names:
                if comp_name in p:
                    para_names.append(p)
            spectral_models.append(SpectralModelConverted(comp.shape, para_names))
        assert len(spectral_models) <= 1, "Only one spec component supported yet"
        if len(spectral_models) > 0:
            self._spectral_model = spectral_models[0]

    def _convert_spatial_model(self) -> None:
        """
        Convert the spatial model of the source if it is an extended one
        """
        raise NotImplementedError("currently fails")
        # need to adapt same style as spectral
        self._spatial_model = SpatialModelConverted(
            self._source.spatial_shape, frame=self._frame
        )

    def _convert_temporal_model(self) -> None:
        """
        Convert the temporal evolution of the source if available
        """
        raise NotImplementedError("currently fails")
        # need to adapt same style as spectral
        self._temporal_model = TemporalModelConverted(self._source.temporal_shape)

    def _create_skymodel(self) -> None:
        """
        Create the skymodel instance out of the individual components
        """
        self._skymodel = SkyModel(
            spectral_model=self._spectral_model,
            spatial_model=self._spatial_model,
            temporal_model=self._temporal_model,
        )
        self._create_parameter_dict()

    def _update_parameter(self, name, val) -> None:
        """
        Update the skymodel parameters during the sampling process
        """
        # update the parameter dict for this skymodel
        self._parameter_dict[name]["value"] = val
        # update the skymodel itself
        self._skymodel.parameters[name].update_from_dict(self._parameter_dict[name])

    @property
    def skymodel(self):
        """
        Returns the Gammapy skymodel for this source
        """
        return self._skymodel

    @property
    def astromdodels_source(self):
        """
        Returns the original astromodel source
        """
        return self._source


class SpectralModelConverted(SpectralModel):
    """
    Class for converting a spectral astromodel function into
    an gammapy SpectralModel
    """

    def __init__(self, function: Function, para_names: list) -> None:
        """
        :param function: the spectral function, must be an astromodels Function
        :type function: astromodels.functions.function.Function
        :param para_names: list of the parameter names for this component
        :type para_names: list[str]
        """
        assert issubclass(
            type(function), Function
        ), "function must be astromodels function"
        self._astromodel_function = function
        self._source_name = self._astromodel_function.name
        self._para_names = para_names
        self._setup_parameters()

    def _setup_parameters(self):
        """
        Setup the parameters by creating gammapy Parameters and setting
        them as attributes to this class
        """
        paras = []
        # needed later for correctly evaluating the function
        self._mapping = {}
        for k, v in self._astromodel_function.parameters.items():
            vmin = np.nan
            vmax = np.nan
            if v.min_value is not None:
                vmin = v.min_value
            if v.max_value is not None:
                vmax = v.max_value
            # find the correct name
            i = 0
            name = None
            while True and i < len(self._para_names):
                splitted = self._para_names[i].split(".")
                if k == splitted[-1]:
                    name = self._para_names[i]
                    break
                i += 1
            if name is None:
                raise ValueError(f"Parameter name {k} not found in {self._para_names}")
            self._mapping[name] = k
            paras.append(
                Parameter(
                    name=name,
                    value=v.value,
                    unit=v.unit,
                    min=vmin,
                    max=vmax,
                    frozen=not bool(v.free),
                )
            )
            setattr(self, name, paras[-1])
        self.default_parameters = Parameters(paras)

    # todo check return type - likely np.ndarray
    def evaluate(self, *args, **kwargs):
        """
        Evaluates the astromodels function instead of a gammapy one
        """

        kwargs_new = {}
        for k in kwargs.keys():
            if k in self._mapping.keys():
                kwargs_new[self._mapping[k]] = kwargs[k]
            else:
                kwargs_new[k] = kwargs[k]
        return self._astromodel_function.evaluate(args[0], **kwargs_new)


class SpatialModelConverted(SpatialModel):
    """
    Class for converting a spatial astromodels function into
    an gammapy SpatialModel
    """

    def __init__(self, function: Function, frame: str = None) -> None:
        """
        :param function: astromodel function describing the morphology
        :type function: astromodels.functions.function.Function
        :param frame: reference frame of the geometry, defaults to ICRS
        :type frame: str

        """
        log.debug("type of spatial function: " + str(type(function)))
        assert issubclass(
            type(function), Function
        ), "function must be astromodels function"
        self._astromodel_function = function
        if frame is None:
            log.warning("No frame passed (not implemented yet) will use ICRS!")
            frame = "icrs"
        self._frame = frame
        self._setup_parameters()

    def _setup_parameters(self) -> None:
        """
        Sets up all the parameters, same as spectral one
        """
        raise NotImplementedError("Need to adapt naming!")
        paras = []
        for k, v in self._astromodel_function.parameters.items():
            vmin = np.nan
            vmax = np.nan
            if v.min_value is not None:
                vmin = v.min_value
            if v.max_value is not None:
                vmax = v.max_value

            paras.append(
                Parameter(
                    name=k,
                    value=v.value,
                    unit=v.unit,
                    min=vmin,
                    max=vmax,
                    frozen=not bool(v.free),
                )
            )
            setattr(self, k, paras[-1])
        setattr(self, "frame", self._frame)
        self.default_parameters = Parameters(paras)

    # todo check return type
    def evaluate(self, *paras, **kwargs):
        """
        Evaluates astromodels function instead of gammapy one
        """
        raise NotImplementedError("Not yet correctly implemented")
        return self._astromodel_function.evaluate(*paras, **kwargs)


class TemporalModelConverted(TemporalModel):
    def __init__(self, function: Function) -> None:
        raise NotImplementedError("Check how this is handled in gammapy")
        log.debug("type of temporal function: " + str(type(function)))
        assert issubclass(
            type(function), Function
        ), "function must be astromodels function"
        self._astromodel_function = function
        self._setup_parameters()

    def _setup_parameters(self):
        paras = []
        for k, v in self._astromodel_function.parameters.items():
            vmin = np.nan
            vmax = np.nan
            if v.min_value is not None:
                vmin = v.min_value
            if v.max_value is not None:
                vmax = v.max_value

            paras.append(
                Parameter(
                    name=k,
                    value=v.value,
                    unit=v.unit,
                    min=vmin,
                    max=vmax,
                    frozen=not bool(v.free),
                )
            )
            setattr(self, k, paras[-1])
        self.default_parameters = Parameters(paras)

    def evaluate(self, *paras, **kwargs):
        return self._astromodel_function.evaluate(*paras, **kwargs)
