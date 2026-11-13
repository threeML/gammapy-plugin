import logging
from copy import deepcopy

import astropy.units as u
import numpy as np
from astromodels.core.parameter_transformation import LogarithmicTransformation
from astromodels.core.sky_direction import SkyDirection
from astromodels.functions.function import Function
from gammapy.modeling.models import (
    PointSpatialModel,
    SpatialModel,
)
from gammapy.modeling.parameter import Parameters

from gammapy_plugin.utils.gammapy_parser import _conversion_list

from .parameters import ParameterWrapped

log = logging.getLogger(__name__)


class SpatialModelConverted(SpatialModel):
    """Class for converting a spatial astromodels function into an gammapy
    SpatialModel."""

    tag = ["SpatialModelConverted", "spat_conv"]

    def __init__(
        self,
        function: Function,
        frame: str | None = None,
    ) -> None:
        """
        :param function: astromodel function describing the morphology
        :param frame: reference frame of the geometry, defaults to ICRS

        """
        log.debug("type of spatial function: " + str(type(function)))
        assert issubclass(type(function), Function), (
            "function must be astromodels function"
        )
        self._astromodel_function = function
        # self._source_name = self._astromodel_function.name
        if frame is None:
            log.info("No frame passed - will use ICRS!")
            frame = "icrs"
        self._frame = frame
        self.frame = self._frame
        self._setup_parameters()
        super().__init__()
        for name, am_parameter in self._am_parameter_mapping.items():
            self.parameters[name].bind(am_parameter)

    def _setup_parameters(self):
        """Setup the parameters by creating gammapy Parameters and setting them
        as attributes to this class."""
        paras = []
        # needed later for correctly evaluating the function
        self._mapping = {}
        self._mapping_free = {}
        self._am_parameter_mapping = {}
        for v in self._astromodel_function.parameters.values():
            factor, unit = _conversion_list.get(v.unit, (1.0, v.unit))
            vmin = v.min_value * factor if v.min_value is not None else np.nan
            vmax = v.max_value * factor if v.max_value is not None else np.nan
            self._mapping[v.path] = v.path
            if v.free:
                self._mapping_free[v.path] = v.path
            is_log = v.has_transformation() and isinstance(
                v.transformation, LogarithmicTransformation
            )
            par = ParameterWrapped(
                name=v.path,
                value=v.value * factor,
                unit=unit,
                min=vmin,
                max=vmax,
                frozen=not v.free,
                interp="log" if is_log else "lin",
            )
            paras.append(par)
            self._am_parameter_mapping[v.path] = v
            setattr(self, v.path, par)
        self.default_parameters = Parameters(paras)

    # todo check return type
    def evaluate(self, *args, **kwargs):
        """Evaluates astromodels function instead of gammapy one."""
        return self._astromodel_function(*args)

    @property
    def mapping(self):
        return self._mapping

    @property
    def mapping_free(self):
        return self._mapping_free


class PointSourceModelConverted(PointSpatialModel):
    tag = ["PointSourceModelConverted", "ps_conv"]

    def __init__(self, sky_position: SkyDirection, frame: str):
        assert isinstance(sky_position, SkyDirection), (
            "sky_position must be SkyDirection"
        )
        self._sky_position = sky_position
        self._name = self._sky_position.name
        self._position = self._sky_position.sky_coord.transform_to(frame)
        self._frame = frame
        log.debug(f"PointSpatialMpdel got frame {self._frame}")
        self.frame = self._frame
        self._setup_parameters()
        super().__init__()

    def _setup_parameters(self):
        """Setup the parameters by creating gammapy Parameters and setting them
        as attributes to this class."""
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
        if lon.value > 180:
            lon -= 360 * u.deg
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
        lon_0 = ParameterWrapped(
            name="lon_0", value=lon.value, unit=lon.unit, frozen=not lon_free
        )
        lat_0 = ParameterWrapped(
            name="lat_0", value=lat.value, unit=lat.unit, frozen=not lat_free
        )
        self.lon_0 = lon_0
        self.lat_0 = lat_0
        self.default_parameters = Parameters([lon_0, lat_0])
        log.debug(f"Set parameters to be {lon_0} and {lat_0}")

    @property
    def mapping(self):
        return self._mapping

    @property
    def mapping_free(self):
        return self._mapping_free

    def _rebind_parameters(self):
        """Bind Gammapy parameters to the actual evaluation function."""

        functions = self._sky_position

        # Get parameters from the function actually used
        # by evaluate(), not from an independently copied
        # astromodels model or parameter mapping.

        if not hasattr(self, "_inv_mapping"):
            self._inv_mapping = {v: k for k, v in self._mapping.items()}

        am_parameters = {par.path: par for par in functions.parameters.values()}

        for gp_parameter in self.parameters:
            if isinstance(gp_parameter, ParameterWrapped):
                am_parameter = am_parameters[self._inv_mapping[gp_parameter.name]]
                gp_parameter.bind(am_parameter)

        # Keep the cached parameter list consistent too.
        self._components_parameters = list(am_parameters.values())

    def __deepcopy__(self, memo):
        """Deep-copy model and restore astromodels parameter bindings."""

        # Respect already-copied objects.
        if id(self) in memo:
            return memo[id(self)]

        # Create an uninitialized instance.
        cls = type(self)
        copied = cls.__new__(cls)

        # Register immediately to handle circular references.
        memo[id(self)] = copied

        # Copy all instance attributes using the SAME memo.
        for name, value in self.__dict__.items():
            object.__setattr__(
                copied,
                name,
                deepcopy(value, memo),
            )

        # Restore parameter references to the copied
        # astromodels function actually used by evaluate().
        copied._rebind_parameters()

        return copied
