import logging

import numpy as np
from astromodels.functions.function import Function
from gammapy.modeling.models import (
    SpectralModel,
)
from gammapy.modeling.parameter import Parameter, Parameters

from gammapy_plugin.utils.gammapy_parser import _conversion_list

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
                _conversion_list[temp_x][1]
                if temp_x in _conversion_list.keys()
                else temp_x
            )
            self._y_unit = (
                _conversion_list[temp_y][1]
                if temp_y in _conversion_list.keys()
                else temp_y
            )

            self._components_parameters = list(
                self._astromodel_function.parameters.values()
            )
        elif isinstance(function, list):
            for f in function:
                assert isinstance(
                    f, Function
                ), f"{f.name} is not a function but {type(f)}"
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
                            if f.x_unit in _conversion_list.keys()
                            else f.x_unit
                        )
                    elif x_unit is not None and f.x_unit is not None:
                        assert x_unit.is_equivalent(
                            f.x_unit
                        ), "Component x_unit not matching"
                        # TODO: maybe transform also possible need to check
                    else:
                        raise ValueError(f"Your Component {f.name} has no x_unit")
                else:
                    x_unit = (
                        _conversion_list[f._requested_x_unit][1]
                        if f._requested_x_unit in _conversion_list.keys()
                        else f._requested_x_unit
                    )
                if not hasattr(f, "_requested_y_unit"):
                    if y_unit is None and f.y_unit is not None:
                        y_unit = (
                            _conversion_list[f.y_unit][1]
                            if f.y_unit in _conversion_list.keys()
                            else f.y_unit
                        )
                    elif y_unit is not None and f.y_unit is not None:
                        assert y_unit.is_equivalent(
                            f.y_unit
                        ), "Component y_unit not matching"
                        # TODO: maybe transform also possible need to check
                    else:
                        raise ValueError(f"Your Component {f.name} has no y_unit")
                else:
                    y_unit = (
                        _conversion_list[f._requested_y_unit][1]
                        if f._requested_y_unit in _conversion_list.keys()
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

    def _setup_parameters(self):
        """Setup the parameters by creating gammapy Parameters and setting them
        as attributes to this class."""
        paras = []
        # needed later for correctly evaluating the function
        self._mapping = {}
        self._mapping_free = {}
        parameter_dict = self._components_parameters
        for v in parameter_dict:
            if v.unit in _conversion_list.keys():
                unit = _conversion_list[v.unit][1]
                factor = _conversion_list[v.unit][0]
            else:
                unit = v.unit
                factor = 1.0
            vmin = np.nan
            vmax = np.nan
            if v.min_value is not None:
                vmin = v.min_value
            if v.max_value is not None:
                vmax = v.max_value
            self._mapping[v.path] = v.path
            if v.free:
                self._mapping_free[v.path] = v.path
            paras.append(
                Parameter(
                    name=v.path,
                    value=v.value * factor,
                    unit=unit,
                    min=vmin * factor,
                    max=vmax * factor,
                    frozen=not bool(v.free),
                )
            )
            setattr(self, v.path, paras[-1])
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
