from typing import Optional
from threeML.io.logging import setup_logger
from gammapy.modeling.models import SkyModel
from astromodels.core.model import Model
from astromodels.sources import PointSource, ExtendedSource, Source
from gammapy_plugin.gammapy_source import parameter_to_gammapy_dict
from gammapy_plugin.models import (
    SpectralModelConverted,
    SpatialModelConverted,
    TemporalModelConverted,
    PointSourceModelConverted,
)

log = setup_logger(__name__)


class AstromodelConverter:
    """
    Class for analyizing an astromodel model and converting all the individual
    sources such that it can be used with gammapy

    Every Source in the Model will get its own Gammapy skymodel.
    The evaluation happens via the astromodel definition.


    :param model: the astromodel model
    :param frame: geometry frame for gammapy, defaults to ICRS

    """

    def __init__(self, model: Model, frame: Optional[str] = None) -> None:
        assert isinstance(model, Model), "Needs an astromodels Model"
        self._astromodel_model = model
        self._frame = frame

        self._converted_sources = {}
        self._gammapy_models = []
        self._convert_extendend_sources()
        self._convert_point_sources()
        self._create_gammapy_models_list()

    def _convert_extendend_sources(self) -> None:
        """
        Converts an extended source into a gammapy skymodel
        """
        for (
            source_name,
            source_instance,
        ) in self._astromodel_model.extended_sources.items():
            self._converted_sources[source_name] = SourceConverter(
                source_instance, frame=self._frame, converter=self
            )

    def _convert_point_sources(self) -> None:
        """
        Converts point sources into individual skymodels
        """
        for (
            source_name,
            source_instance,
        ) in self._astromodel_model.point_sources.items():
            self._converted_sources[source_name] = SourceConverter(
                source_instance, converter=self
            )

    def _create_gammapy_models_list(self) -> None:
        """
        Creates a list with all SkyModels fro the converted sources
        """
        for name, source in self._converted_sources.items():
            self._gammapy_models.append(source.skymodel)

    def _update_parameters(self) -> None:
        """
        Update all the parameters in the SkyModels with the values from
        the astromodels model
        """
        # TODO this is a very stupid way of checking if a parameter
        # belongs to a source
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
    """
    Takes a astromodels source and converts it to a SkyModel
    :param source: The astromodels source
    :param converter: The used AstromodelConverter instance
    :param kwargs: frame
    """

    def __init__(
        self, source: Source, converter: AstromodelConverter, **kwargs
    ) -> None:
        self._frame = kwargs.get("frame", "icrs")
        self._converter = converter
        self._spatial_model = None
        self._spectral_model = None
        self._temporal_model = None
        self._gammapy_model = None
        self._parameter_dict = None
        self._source = source
        self._all_para_names = list(self._source.parameters.keys())

        self._convert_spatial_model()
        self._convert_spectral_model()

        self._create_skymodel()

    def _create_parameter_dict(self) -> None:
        """
        Initially creates the parameter dict for that source
        the SkyModel is updated using the update_from_dict function
        """

        # TODO thats likely inefficient
        # TODO use mapping from Models
        if self._parameter_dict is None:
            self._parameter_dict = {}

        for name in self._skymodel.parameters.names:
            self._parameter_dict[name] = {}
            if name in self._converter._astromodel_model.parameters.keys():
                astromodel_para = self._converter._astromodel_model.parameters[name]
            elif (
                self._source.name + "." + name
                in self._converter._astromodel_model.parameters.keys()
            ):
                astromodel_para = self._converter._astromodel_model.parameters[
                    self._source.name + "." + name
                ]
            elif (
                self._source.name + "." + self._source.name + "." + name
                in self._converter._astromodel_model.parameters.keys()
            ):
                astromodel_para = self._converter._astromodel_model.parameters[
                    self._source.name + "." + self._source.name + "." + name
                ]
            else:
                log.error(f"The skymodel parameter {name} is not known")
                log.error(
                    f"These are the astromodel paras {self._converter._astromodel_model.parameters.keys()}"
                )
                raise NotImplementedError
            self._parameter_dict[name] = parameter_to_gammapy_dict(astromodel_para)

    def _update_parameter(self, name, val) -> None:
        """
        Update the skymodel parameters during the sampling process
        using the parameter dict
        """
        # update the parameter dict for this skymodel
        # update the skymodel itself
        self._parameter_dict[name]["value"] = val
        self._skymodel.parameters[name].update_from_dict(self._parameter_dict[name])

    def _convert_spectral_model(self) -> None:
        """
        Convert the spectral model of the source
        """
        log.warning("Multiple spectral components well be simply added!")
        spectral_models = []
        for comp_name, comp in self._source.components.items():
            para_names = []
            for p in self._all_para_names:
                if comp_name in p:
                    para_names.append(p)
            spectral_models.append(SpectralModelConverted(comp.shape, para_names))
        self._spectral_model = spectral_models[0]
        if len(spectral_models) > 1:
            for spectral_model in spectral_models:
                self._spectral_model += spectral_model

    def _convert_spatial_model(self) -> None:
        """
        Convert the spatial model of the source
        """
        if isinstance(self._source, ExtendedSource):
            comp_name = self._source.spatial_shape.name
            ps = False
        elif isinstance(self._source, PointSource):
            comp_name = "position"
            position = self._source.position
            ps = True
        para_names = []
        for p in self._all_para_names:
            if comp_name in p:
                para_names.append(p)
        if ps:
            self._spatial_model = PointSourceModelConverted(
                sky_position=position,
                frame=self._frame,
                para_names=para_names,
            )
        else:
            self._spatial_model = SpatialModelConverted(
                self._source.spatial_shape, para_names, frame=self._frame
            )

    def _convert_temporal_model(self) -> None:
        """
        Convert the temporal evolution of the source if available
        """
        raise NotImplementedError("Not yet implemented")
        # need to adapt same style as spectral
        self._temporal_model = TemporalModelConverted(self._source.temporal_shape)

    def _create_skymodel(self) -> None:
        """
        Create the skymodel instance out of the individual components
        """
        if self._gammapy_model is None:
            self._skymodel = SkyModel(
                name=self._source.name,
                spectral_model=self._spectral_model,
                spatial_model=self._spatial_model,
                temporal_model=self._temporal_model,
            )
        else:
            self._skymodel = self._gammapy_model
        self._create_parameter_dict()

    @property
    def skymodel(self) -> SkyModel:
        """
        Returns the Gammapy skymodel for this source
        """
        return self._skymodel

    @property
    def astromdodels_source(self) -> Source:
        """
        Returns the original astromodel source
        """
        return self._source
