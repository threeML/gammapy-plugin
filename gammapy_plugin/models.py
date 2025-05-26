import astropy.units as u
import numpy as np
from astromodels.core.parameter_transformation import LogarithmicTransformation
from astromodels.core.sky_direction import SkyDirection
from astromodels.functions.function import Function
from gammapy.modeling.covariance import Covariance
from gammapy.modeling.models import (
    PointSpatialModel,
    SpatialModel,
    SpectralModel,
    TemporalModel,
)
from gammapy.modeling.parameter import Parameter, Parameters
from threeML.io.logging import setup_logger

log = setup_logger(__name__)


class SpectralModelConverted(SpectralModel):
    """
    Class for converting a spectral astromodel function into
    an gammapy SpectralModel
    """

    tag = ["SpectralModelConverted", "spec_conv"]

    def __init__(self, function: Function) -> None:
        """
        :param function: the spectral function, must be an astromodels Function
        """
        assert issubclass(
            type(function), Function
        ), "function must be astromodels function"
        self._astromodel_function = function
        # self._source_name = self._astromodel_function.name
        self._setup_parameters()
        self._x_unit = self._astromodel_function.x_unit
        self._y_unit = self._astromodel_function.y_unit
        log.debug(f"These are the units: {self._x_unit}, {self._y_unit}")
        self._integral_unit = self._y_unit * self._x_unit
        super().__init__()

    def _setup_parameters(self):
        """
        Setup the parameters by creating gammapy Parameters and setting
        them as attributes to this class
        """
        paras = []
        # needed later for correctly evaluating the function
        self._mapping = {}
        self._mapping_free = {}
        for k, v in self._astromodel_function.parameters.items():
            vmin = np.nan
            vmax = np.nan
            if v.min_value is not None:
                vmin = v.min_value
            if v.max_value is not None:
                vmax = v.max_value
            self._mapping[v.path] = v.name
            if v.free:
                self._mapping_free[v.path] = v.name
            # interp = "linear"
            # if isinstance(v.transformation, LogarithmicTransformation):
            #    interp = "log"
            #    log.debug(f"Set interpolatin of {v.name} to logarithmic")
            paras.append(
                Parameter(
                    name=v.name,
                    value=v.value,
                    unit=v.unit,
                    min=vmin,
                    max=vmax,
                    frozen=not bool(v.free),
                    # interp=interp,
                    # scale_transform=interp,
                )
            )
            setattr(self, v.name, paras[-1])
        self.default_parameters = Parameters(paras)

    def evaluate(self, *args, **kwargs):
        """
        Evaluates the astromodels function instead of a gammapy one
        """
        # TODO check the kwarg stuff
        kwargs_new = {}
        for k in kwargs.keys():
            if k in self._mapping.keys():
                kwargs_new[self._mapping[k]] = kwargs[k]
            else:
                kwargs_new[k] = kwargs[k]
        return self._astromodel_function.evaluate(*args, **kwargs_new)

    @property
    def mapping(self):
        return self._mapping

    @property
    def mapping_free(self):
        return self._mapping_free


class PointSourceModelConverted(PointSpatialModel):
    tag = ["PointSourceModelConverted", "ps_conv"]

    def __init__(self, sky_position: SkyDirection, frame: str):
        assert isinstance(
            sky_position, SkyDirection
        ), "sky_position must be SkyDirection"
        self._sky_position = sky_position
        self._name = self._sky_position.name
        self._position = self._sky_position.sky_coord.transform_to(frame)
        self._frame = frame
        log.debug(f"PointSpatialMpdel got frame {self._frame}")
        setattr(self, "frame", self._frame)
        self._setup_parameters()
        super().__init__()

    def _setup_parameters(self):
        """
        Setup the parameters by creating gammapy Parameters and setting
        them as attributes to this class
        """
        self._mapping = {}
        self._mapping_free = {}
        if self._frame == "galactic":
            lon = self._position.l
            lat = self._position.b
        elif self._frame == "icrs":
            lon = self._position.ra
            lat = self._position.dec
        else:
            raise NotImplementedError("Only galactic and icrs currently available")
        # TODO check if this is necessary
        if lon.value > 180:
            lon -= 360 * u.deg
        elif lon.value <= -180:
            lon += 360 * u.deg
        for k, v in self._sky_position.parameters.items():
            if k == "ra" or k == "l":
                para_name = "lon_0"
                lon_free = v.free
            elif k == "dec" or k == "b":
                para_name = "lat_0"
                lat_free = v.free
            self._mapping[v.path] = para_name
            if v.free:
                self._mapping_free[v.path] = para_name
        lon_0 = Parameter(
            name="lon_0", value=lon.value, unit=lon.unit, frozen=not lon_free
        )
        lat_0 = Parameter(
            name="lat_0", value=lat.value, unit=lat.unit, frozen=not lat_free
        )
        setattr(self, "lon_0", lon_0)
        setattr(self, "lat_0", lat_0)
        self.default_parameters = Parameters([lon_0, lat_0])
        log.debug(f"Set parameters to be {lon_0} and {lat_0}")

    @property
    def mapping(self):
        return self._mapping

    @property
    def mapping_free(self):
        return self._mapping_free


class SpatialModelConverted(SpatialModel):
    """
    Class for converting a spatial astromodels function into
    an gammapy SpatialModel
    """

    tag = ["SpatialModelConverted", "spat_conv"]

    def __init__(
        self,
        function: Function,
        frame: str = None,
    ) -> None:
        """
        :param function: astromodel function describing the morphology
        :param frame: reference frame of the geometry, defaults to ICRS

        """
        log.debug("type of spatial function: " + str(type(function)))
        assert issubclass(
            type(function), Function
        ), "function must be astromodels function"
        self._astromodel_function = function
        # self._source_name = self._astromodel_function.name
        if frame is None:
            log.warning("No frame passed - will use ICRS!")
            frame = "icrs"
        self._frame = frame
        setattr(self, "frame", self._frame)
        self._setup_parameters()
        super().__init__()

    def _setup_parameters(self):
        """
        Setup the parameters by creating gammapy Parameters and setting
        them as attributes to this class
        """
        paras = []
        # needed later for correctly evaluating the function
        self._mapping = {}
        self._mapping_free = {}
        for k, v in self._astromodel_function.parameters.items():
            vmin = np.nan
            vmax = np.nan
            if v.min_value is not None:
                vmin = v.min_value
            if v.max_value is not None:
                vmax = v.max_value
            self._mapping[v.path] = v.name
            if v.free:
                self._mapping_free[v.path] = v.name
            interp = "linear"
            if isinstance(v.transformation, LogarithmicTransformation):
                interp = "log"
            paras.append(
                Parameter(
                    name=v.name,
                    value=v.value,
                    unit=v.unit,
                    min=vmin,
                    max=vmax,
                    frozen=not bool(v.free),
                    interp=interp,
                )
            )
            setattr(self, v.name, paras[-1])
        self.default_parameters = Parameters(paras)

    # todo check return type
    def evaluate(self, *args, **kwargs):
        """
        Evaluates astromodels function instead of gammapy one
        """
        kwargs_new = {}
        for k in kwargs.keys():
            if k in self._mapping.keys():
                kwargs_new[self._mapping[k]] = kwargs[k]
            else:
                kwargs_new[k] = kwargs[k]
        return self._astromodel_function.evaluate(*args, **kwargs_new)

    @property
    def mapping(self):
        return self._mapping

    @property
    def mapping_free(self):
        return self._mapping_free


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

    @property
    def mapping(self):
        return self._mapping
