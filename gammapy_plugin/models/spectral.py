import logging
from copy import deepcopy

import numpy as np
from astromodels.core.parameter_transformation import LogarithmicTransformation
from astromodels.functions.function import Function
from gammapy.modeling.models import (
    SpectralModel,
)
from gammapy.modeling.parameter import Parameters

from gammapy_plugin.utils.gammapy_parser import _conversion_list

from .parameters import ParameterWrapped

log = logging.getLogger(__name__)


class SpectralModelConverted(SpectralModel):
    """Class for converting a spectral astromodel function into an gammapy
    SpectralModel."""

    tag = ["SpectralModelConverted", "spec_conv"]

    def __init__(self, function: Function | list, **kwargs) -> None:

        self._components_parameters = None
        if isinstance(function, Function):
            self._astromodel_function = function
            self._components = None
            if hasattr(self._astromodel_function, "_requested_x_unit"):
                # this is a composite function
                temp_x = self._astromodel_function._requested_x_unit
                temp_y = self._astromodel_function._requested_y_unit
            else:
                temp_x = self._astromodel_function.x_unit
                temp_y = self._astromodel_function.y_unit

            self._x_unit = (
                _conversion_list[temp_x][1] if temp_x in _conversion_list else temp_x
            )
            self._y_unit = (
                _conversion_list[temp_y][1] if temp_y in _conversion_list else temp_y
            )

            self._components_parameters = list(
                self._astromodel_function.parameters.values()
            )
        elif isinstance(function, list):
            for f in function:
                assert isinstance(f, Function), (
                    f"{f.name} is not a function but {type(f)}"
                )
            self._astromodel_function = function
            self._components = len(function)
            x_unit = None
            y_unit = None

            self._components_parameters = []
            for f in self._astromodel_function:
                if not hasattr(f, "_requested_x_unit"):
                    if x_unit is None and f.x_unit is not None:
                        x_unit = (
                            _conversion_list[f.x_unit][1]
                            if f.x_unit in _conversion_list
                            else f.x_unit
                        )
                    elif x_unit is not None and f.x_unit is not None:
                        assert x_unit.is_equivalent(f.x_unit), (
                            "Component x_unit not matching"
                        )
                        # TODO: maybe transform also possible need to check
                    else:
                        raise ValueError(f"Your Component {f.name} has no x_unit")
                else:
                    x_unit = (
                        _conversion_list[f._requested_x_unit][1]
                        if f._requested_x_unit in _conversion_list
                        else f._requested_x_unit
                    )
                if not hasattr(f, "_requested_y_unit"):
                    if y_unit is None and f.y_unit is not None:
                        y_unit = (
                            _conversion_list[f.y_unit][1]
                            if f.y_unit in _conversion_list
                            else f.y_unit
                        )
                    elif y_unit is not None and f.y_unit is not None:
                        assert y_unit.is_equivalent(f.y_unit), (
                            "Component y_unit not matching"
                        )
                        # TODO: maybe transform also possible need to check
                    else:
                        raise ValueError(f"Your Component {f.name} has no y_unit")
                else:
                    y_unit = (
                        _conversion_list[f._requested_y_unit][1]
                        if f._requested_y_unit in _conversion_list
                        else f._requested_y_unit
                    )

                for p in f.parameters.values():
                    self._components_parameters.append(p)
            self._x_unit = x_unit
            self._y_unit = y_unit

        else:
            raise NotImplementedError(
                "Can only convert astromodels Function or list of Functions"
            )

        if self._x_unit is None or self._y_unit is None:
            raise ValueError("You need to specify units for your spectral component")
        log.debug(f"These are the units: {self._x_unit}, {self._y_unit}")

        self._setup_parameters()

        self._integral_unit = self._y_unit * self._x_unit
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
        for v in self._components_parameters:
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

    def evaluate(self, energy, **kwargs):
        """Evaluates the astromodels function instead of a gammapy one."""
        shape = None
        if len(energy.shape) > 1:
            shape = energy.shape
            energy = energy.flatten()
        if self._components is not None:
            vals = []
            if shape is None:
                for i in range(self._components):
                    val = self._astromodel_function[i](energy)
                    vals.append(val)
            else:
                for i in range(self._components):
                    vals.append(self._astromodel_function[i](energy).reshape(shape))

            return sum(vals)

        else:
            if shape is None:
                return self._astromodel_function(energy)
            else:
                return self._astromodel_function(energy).reshape(shape)

    @property
    def mapping(self):
        return self._mapping

    @property
    def mapping_free(self):
        return self._mapping_free

    def _rebind_parameters(self):
        """Bind Gammapy parameters to the actual evaluation function."""

        functions = self._astromodel_function

        if not isinstance(functions, list):
            functions = [functions]

        # Get parameters from the function actually used
        # by evaluate(), not from an independently copied
        # astromodels model or parameter mapping.
        am_parameters = {
            par.path: par
            for function in functions
            for par in function.parameters.values()
        }

        for gp_parameter in self.parameters:
            if isinstance(gp_parameter, ParameterWrapped):
                am_parameter = am_parameters[gp_parameter.name]
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
