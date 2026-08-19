import logging

from astromodels.core.model import Model
from astromodels.sources import ExtendedSource, PointSource, Source
from gammapy.modeling.models import SkyModel

from gammapy_plugin.models import (
    PointSourceModelConverted,
    SpatialModelConverted,
    SpectralModelConverted,
)
from gammapy_plugin.utils.gammapy_parser import (
    parameter_to_gammapy_dict,
)

__all__ = ["AstromodelConverter", "SourceConverter"]

log = logging.getLogger(__name__)


class AstromodelConverter:
    """Class for analyizing an astromodel model and converting all the
    individual sources such that it can be used with gammapy.

    Every Source in the Model will get its own Gammapy skymodel. The
    evaluation happens via the astromodel definition.
    """

    def __init__(
        self, model: Model, frame: str | None = None, convert_ps: bool = True
    ) -> None:
        """
        :param model: the astromodel model
        :type model: astromodels.core.model.Model
        :param frame: geometry frame for gammapy, defaults to ICRS, optional
        :type frame: str
        :param convert_ps: convert PointSources to a PointSpatialModel, optional
        :type convert_ps: bool
        """
        assert isinstance(
            model, Model
        ), "AstromodelConverter needs an astromodels Model"
        self._astromodel_model = model
        if frame is not None:
            self._frame = frame
        else:
            log.info("No frame passed - will use ICRS")
            self._frame = "icrs"
        self._convert_ps = convert_ps
        self._converted_sources = {}
        self._gammapy_models = []
        self._convert_extendend_sources()
        self._convert_point_sources()
        self._create_gammapy_models_list()

    def _convert_extendend_sources(self) -> None:
        """Converts all extended sources into individual gammapy skymodels."""
        for (
            source_name,
            source_instance,
        ) in self._astromodel_model.extended_sources.items():
            self._converted_sources[source_name] = SourceConverter(
                source_instance, frame=self._frame, converter=self
            )

    def _convert_point_sources(self) -> None:
        """Converts all point sources into individual skymodels."""
        for (
            source_name,
            source_instance,
        ) in self._astromodel_model.point_sources.items():
            log.debug(f"Converting pointsource {source_name}")
            self._converted_sources[source_name] = SourceConverter(
                source_instance, converter=self, convert_spatial=self._convert_ps
            )

    def _create_gammapy_models_list(self) -> None:
        """Creates a list with all SkyModels fro the converted sources."""
        for name, source in self._converted_sources.items():
            self._gammapy_models.append(source.skymodel)

    def _update_parameters(self) -> None:
        """Update all the parameters in the SkyModels with the values from the
        astromodels model using the underlying SourceConverters."""
        for name, source in self._converted_sources.items():
            source._update_parameters()

    def update(self) -> None:
        """Update all parameters of all SkyModels witht the current values from the
        astromodels model. Public method for _update_parameters"""
        self._update_parameters()

    @property
    def gammapy_models(self) -> list[SkyModel]:
        """
        Returns all the gammapy skymodels for that model.

        :return: list with all the SkyModels
        :rtype: list[SkyModel]
        """
        return self._gammapy_models

    @property
    def model(self) -> Model:
        """
        Return the corresponding astromodels model

        :return: the astromodels model
        :rtype: astromodels.core.model.Model
        """

        return self._astromodel_model

    @property
    def converted_sources(self) -> dict:
        """
        Returns dictionary with all the converted sources.

        :return: the converted sources dict
        :rtype: dict
        """
        return self._converted_sources


class SourceConverter:
    """Takes a astromodels source and converts it to a SkyModel"""

    def __init__(
        self, source: Source, converter: AstromodelConverter, **kwargs
    ) -> None:
        """
        Initialize the SourceConverter.
        :param source: The astromodels source
        :type source: astromodels.core.model.Model
        :param converter: The used AstromodelConverter instance
        :type converter: AstromodelConverter
        :param kwargs: Additional keyword arguments, optional
        :type kwrags: dict
        """

        self._conv_spatial = kwargs.get("convert_spatial", True)
        self._converter = converter
        self._frame = self._converter._frame
        self._source = source

        # init everything to None
        self._spatial_model = None
        self._spectral_model = None
        self._temporal_model = None
        self._gammapy_model = None
        self._parameter_dict = None

        # start with converting the spatial model
        if self._conv_spatial:
            log.debug("Converting the spatial Model")
            self._convert_spatial_model()

        log.debug("Converting the spectral Model")
        # converting the spectral model
        self._convert_spectral_model()

        self._create_skymodel()

    def _gather_mappings(self) -> None:
        """
        Gather all the parameter mappings from astromodels to gammapy for the individual
        spectral and spatial models
        """
        self._mapping = None
        self._mapping_free = None
        for comp in [self._spectral_model, self._spatial_model]:
            if comp is not None:
                if self._mapping is None:
                    self._mapping = comp.mapping
                else:
                    self._mapping.update(comp.mapping)
                if self._mapping_free is None:
                    self._mapping_free = comp.mapping_free
                else:
                    self._mapping_free.update(comp.mapping_free)

    def _update_parameters(self) -> None:
        """Update the skymodel parameters during the sampling process using the
        parameter dict."""
        # update the parameter dict for this skymodel
        for k, v in self._mapping_free.items():
            self.skymodel.parameters[v].update_from_dict(
                parameter_to_gammapy_dict(self._converter.model[k])
            )

    def _convert_spectral_model(self) -> None:
        """Convert the spectral model of the source."""
        if len(self._source.components.keys()) > 1:
            self._spectral_model = SpectralModelConverted(
                [x.shape for x in self._source.components.values()],
            )

        else:
            self._spectral_model = SpectralModelConverted(
                self._source.components[list(self._source.components.keys())[0]].shape,
            )

    def _convert_spatial_model(self) -> None:
        """Convert the spatial model of the source."""
        if isinstance(self._source, ExtendedSource):
            ps = False
        elif isinstance(self._source, PointSource):
            position = self._source.position
            ps = True
        if ps:
            self._spatial_model = PointSourceModelConverted(
                sky_position=position,
                frame=self._frame,
            )
        else:
            self._spatial_model = SpatialModelConverted(
                self._source.spatial_shape, frame=self._frame
            )

    def _create_skymodel(self) -> None:
        """Create the skymodel instance out of the individual components."""
        self._skymodel = SkyModel(
            name=self._source.name,
            spectral_model=self._spectral_model,
            spatial_model=self._spatial_model,
            temporal_model=self._temporal_model,
        )
        self._gather_mappings()

    def update(self) -> None:
        """
        This invokes the updating of all parameters of the converted model in the
        gammapy SkyModels
        """
        self._update_parameters()

    @property
    def skymodel(self) -> SkyModel:
        """Returns the Gammapy skymodel for this source.

        :returns: the created gammapy SkyModel
        :rtype: gammapy.modeling.models.SkyModel

        """
        return self._skymodel

    @property
    def astromodels_source(self) -> Source:
        """Returns the original astromodel source.

        :returns: the converted astromodels Source
        :rtype: astromodels.sources.Source
        """

        return self._source


class GammapyConverter:  # pragma: no cover
    """Takes a gammapy SkyModel and transforms it into a astromodels Model."""

    def __init__(self, model: SkyModel) -> None:
        assert isinstance(model, SkyModel), "You must provide a gammapy.model.SkyModel"
        self._sky_model = model
        raise NotImplementedError("Sorry this is not yet a feature :(")

    def _parse_sources(self):
        spectral_model = self._sky_model.spectral_model  # noqa: F841
        spatial_model = self._sky_model.spatial_model  # noqa: F841

    def _create_astromodel(self):
        # load the sources
        # load the
        sources = []
        self._astromodel = Model(*sources)
