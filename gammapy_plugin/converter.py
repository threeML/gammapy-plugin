import logging
from typing import Optional

from astromodels.core.model import Model
from astromodels.sources import ExtendedSource, PointSource, Source
from gammapy.modeling.models import SkyModel

from gammapy_plugin.gammapy_source import parameter_to_gammapy_dict
from gammapy_plugin.models import (
    PointSourceModelConverted,
    SpatialModelConverted,
    SpectralModelConverted,
)

log = logging.getLogger(__name__)


class AstromodelConverter:
    """Class for analyizing an astromodel model and converting all the
    individual sources such that it can be used with gammapy.

    Every Source in the Model will get its own Gammapy skymodel. The
    evaluation happens via the astromodel definition.

    :param model: the astromodel model
    :param frame: geometry frame for gammapy, defaults to ICRS
    """

    def __init__(
        self, model: Model, frame: Optional[str] = None, convert_ps: bool = True
    ) -> None:
        assert isinstance(model, Model), "Needs an astromodels Model"
        self._astromodel_model = model
        if frame is not None:
            self._frame = frame
        else:
            log.warning("No frame passed - will use ICRS")
            self._frame = "icrs"
        self._convert_ps = convert_ps
        self._converted_sources = {}
        self._gammapy_models = []
        self._convert_extendend_sources()
        self._convert_point_sources()
        self._create_gammapy_models_list()

    def _convert_extendend_sources(self) -> None:
        """Converts an extended source into a gammapy skymodel."""
        for (
            source_name,
            source_instance,
        ) in self._astromodel_model.extended_sources.items():
            self._converted_sources[source_name] = SourceConverter(
                source_instance, frame=self._frame, converter=self
            )

    def _convert_point_sources(self) -> None:
        """Converts point sources into individual skymodels."""
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
        astromodels model."""
        # TODO this is a very stupid way of checking if a parameter
        # belongs to a source
        for name, source in self._converted_sources.items():
            source._update_parameters()

    @property
    def gammapy_models(self) -> list[SkyModel]:
        """Returns all the gammapy skymodels for that model."""
        return self._gammapy_models

    @property
    def model(self) -> Model:
        return self._astromodel_model


class SourceConverter:
    """Takes a astromodels source and converts it to a SkyModel :param source:

    The astromodels source
    :param converter: The used AstromodelConverter instance
    :param kwargs: frame.
    """

    def __init__(
        self, source: Source, converter: AstromodelConverter, **kwargs
    ) -> None:
        self._conv_spatial = kwargs.get("convert_spatial", True)
        self._converter = converter
        self._frame = self._converter._frame
        self._source = source
        self._spatial_model = None
        self._spectral_model = None
        self._temporal_model = None
        self._gammapy_model = None
        self._parameter_dict = None
        self._spatial_correction = False
        if self._conv_spatial:
            log.debug("Converting the spatial Model")
            self._convert_spatial_model()
        log.debug("Converting the spectral Model")
        self._convert_spectral_model()

        self._create_skymodel()

    def _gather_mappings(self):
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
        log.warning("Only Single Spectral Component Models currently supported")
        self._spectral_model = SpectralModelConverted(
            self._source.spectrum._get_children()[0].shape,
            spatial_correction=self._spatial_correction,
        )

    def _convert_spatial_model(self) -> None:
        """Convert the spatial model of the source."""
        if isinstance(self._source, ExtendedSource):
            ps = False
            self._spatial_correction = True
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

    @property
    def skymodel(self) -> SkyModel:
        """Returns the Gammapy skymodel for this source."""
        return self._skymodel

    @property
    def astromodels_source(self) -> Source:
        """Returns the original astromodel source."""
        return self._source
