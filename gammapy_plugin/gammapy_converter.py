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
    Class for analyizing a astromodel model and converting all the individual
    sources
    """

    # TODO need to find way to connect temporal evolution
    def __init__(self, model: Model, frame: str = None) -> None:
        assert isinstance(model, Model), "Needs an astromodels Model"

        self._astromodel_model = model
        self._frame = frame

        self._converted_sources = {}
        self._gammapy_models = []
        self._convert_extendend_sources()
        self._convert_point_sources()
        self._create_gammapy_models()

    def _convert_extendend_sources(self):
        for (
            source_name,
            source_instance,
        ) in self._astromodel_model.extended_sources.items():
            self._converted_sources[source_name] = SourceConverter(
                source_instance, frame=self._frame, converter=self
            )

    def _convert_point_sources(self):
        for (
            source_name,
            source_instance,
        ) in self._astromodel_model.point_sources.items():
            self._converted_sources[source_name] = SourceConverter(
                source_instance, converter=self
            )

    def _create_gammapy_models(self):
        for source in self._converted_sources:
            self._gammapy_models.append(source.skymodel)

    @property
    def gammapy_models(self) -> list[SkyModel]:
        return self._gammapy_models


class SourceConverter:
    def __init__(self, source: Source, **kwargs) -> None:
        self._frame = kwargs.get("frame", None)
        self._spatial_model = None
        self._spectral_model = None
        self._temporal_model = None

        self._source = source
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

    def _convert_spectral_model(self):
        spectral_models = []
        for comp_name, comp in self._source.components.items():
            spectral_models.append(SpectralModelConverted(comp.shape))
        assert len(spectral_models) <= 1, "Only one spec component supported yet"
        if len(spectral_models) > 0:
            self._spectral_model = spectral_models[0]

    def _convert_spatial_model(self):
        self._spatial_model = SpatialModelConverted(
            self._source.spatial_shape, frame=self._frame
        )

    def _convert_temporal_model(self):
        self._temporal_model = TemporalModelConverted(self._source.temporal_shape)

    def _create_skymodel(self):
        self._skymodel = SkyModel(
            spectral_model=self._spectral_model,
            spatial_model=self._spatial_model,
            temporal_model=self._temporal_model,
        )

    @property
    def skymodel(self):
        return self._skymodel

    @property
    def astromdodels_source(self):
        return self._source


class SpectralModelConverted(SpectralModel):
    def __init__(self, function: Function) -> None:
        log.debug("type of spectral function: " + str(type(function)))
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

    def evaluate(self, *args, **kwargs):
        return self._astromodel_function.evaluate(args[0], **kwargs)


class SpatialModelConverted(SpatialModel):
    def __init__(self, function: Function, frame: str = None) -> None:
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
        setattr(self, "frame", self._frame)
        self.default_parameters = Parameters(paras)

    def evaluate(self, *paras, **kwargs):
        return self._astromodel_function.evaluate(*paras, **kwargs)


class TemporalModelConverted(TemporalModel):
    def __init__(self, function: Function) -> None:
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
