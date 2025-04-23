from threeML.io.logging import setup_logger
import numpy as np
import astropy.units as u
from gammapy.modeling.models import (
    SpectralModel,
    SpatialModel,
    TemporalModel,
    PointSpatialModel,
)
from gammapy.modeling.parameter import Parameter, Parameters
from astromodels.functions.function import Function
from astromodels.core.sky_direction import SkyDirection

log = setup_logger(__name__)


class SpectralModelConverted(SpectralModel):
    """
    Class for converting a spectral astromodel function into
    an gammapy SpectralModel
    """

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
        self._integral_unit = self._y_unit * self._x_unit

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
            self._mapping[v.path] = v.name
            paras.append(
                Parameter(
                    name=v.name,
                    value=v.value,
                    unit=v.unit,
                    min=vmin,
                    max=vmax,
                    frozen=not bool(v.free),
                )
            )
            setattr(self, v.name, paras[-1])
        self.default_parameters = Parameters(paras)

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
        return self._astromodel_function.evaluate(*args, **kwargs_new)

    def evaluate_integral(
        self, emin: u.Quantity, emax: u.Quantity, **kwargs
    ) -> u.Quantity:
        """
        Custom integral
        """
        assert len(emin) == len(emax), "Energy edges length do not match"
        vals = np.zeros(len(emin))
        emin = emin.to(self._x_unit).value
        emax = emax.to(self._x_unit).value
        for i in range(len(emin)):
            x = np.linspace(emin[i], emax[i], num=100)
            vals[i] = np.trapezoid(self._astromodel_function.fast_call(x))
        return vals * self._integral_unit


class PointSourceModelConverted(PointSpatialModel):
    def __init__(self, sky_position: SkyDirection, frame: str, para_names: list[str]):
        assert isinstance(
            sky_position, SkyDirection
        ), "sky_position must be SkyDirection"
        self._sky_position = sky_position
        self._name = self._sky_position.name
        self._position = self._sky_position.sky_coord.transform_to(frame)
        self._frame = frame
        setattr(self, "frame", self._frame)
        self._para_names = para_names
        self._setup_parameters()

    def _setup_parameters(self):
        """
        Setup the parameters by creating gammapy Parameters and setting
        them as attributes to this class
        """
        self._mapping = {}
        if self._frame == "galactic":
            lon = self._position.l
            lat = self._position.b
        elif self._frame == "icrs":
            lon = self._position.ra
            lat = self._position.dec
        else:
            raise NotImplementedError("Only galactic and icrs currently available")
        lon_0 = Parameter(name="lon_0", value=lon.value, unit=lon.unit)
        lat_0 = Parameter(name="lat_0", value=lat.value, unit=lat.unit)
        setattr(self, "lon_0", lon_0)
        setattr(self, "lat_0", lat_0)
        self.default_parameters = Parameters([lon_0, lat_0])
        log.debug(f"Set parameters to be {lon_0} and {lat_0}")


class SpatialModelConverted(SpatialModel):
    """
    Class for converting a spatial astromodels function into
    an gammapy SpatialModel
    """

    def __init__(
        self,
        function: Function,
        para_names: list,
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
            # TODO this is fairly inefficient
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
