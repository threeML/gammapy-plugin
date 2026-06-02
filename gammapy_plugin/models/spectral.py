import logging
from typing import Union

import numpy as np
from astromodels.functions.function import Function
from gammapy.modeling.models import (
    SpectralModel,
)
from gammapy.modeling.parameter import Parameter, Parameters

log = logging.getLogger(__name__)


class SpectralModelConverted(SpectralModel):
    """Class for converting a spectral astromodel function into an gammapy
    SpectralModel."""

    tag = ["SpectralModelConverted", "spec_conv"]

    def __init__(self, function: Union[Function, list], **kwargs) -> None:

        self._components_parameters = None
        if isinstance(function, Function):
            self._astromodel_function = function
            self._components = None
            self._x_unit = self._astromodel_function.x_unit
            self._y_unit = self._astromodel_function.y_unit
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
                if x_unit is None and f.x_unit is not None:
                    x_unit = f.x_unit
                elif x_unit is not None and f.x_unit is not None:
                    assert x_unit == f.x_unit, "Component x_unit not matching"
                    # TODO maybe transform also possible need to check
                else:
                    raise ValueError(f"Your Component {f.name} has no x_unit")

                if y_unit is None and f.y_unit is not None:
                    y_unit = f.y_unit
                elif y_unit is not None and f.y_unit is not None:
                    assert y_unit == f.y_unit, "Component y_unit not matching"
                    # TODO maybe transform also possible need to check
                else:  # pragma: no cover
                    raise ValueError(f"Your Component {f.name} has no y_unit")
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
        parameter_dict = (
            self._astromodel_function.parameters.values()
            if self._components_parameters is None
            else self._components_parameters
        )
        for v in parameter_dict:
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
                    value=v.value,
                    unit=v.unit,
                    min=vmin,
                    max=vmax,
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
                    kwargs_mapped = {}
                    for k, v in kwargs.items():
                        if self._astromodel_function[i].path in k:
                            kwargs_mapped[
                                k.split(self._astromodel_function[i].path + ".")[1]
                            ] = v
                    val = self._astromodel_function[i].evaluate(energy, **kwargs_mapped)
                    vals.append(val)
            else:
                for i in range(self._components):
                    kwargs_mapped = {}
                    for k, v in kwargs.items():
                        if self._astromodel_function[i].path in k:
                            kwargs_mapped[
                                k.split(self._astromodel_function[i].path + ".")[1]
                            ] = v

                    vals.append(
                        self._astromodel_function[i]
                        .evaluate(energy, **kwargs_mapped)
                        .reshape(shape)
                    )

            return sum(vals)

        else:
            kwargs_mapped = {}
            for k, v in kwargs.items():
                if self._astromodel_function.path in k:
                    kwargs_mapped[k.split(f"{self._astromodel_function.path}.")[1]] = v
            if shape is None:
                return self._astromodel_function.evaluate(energy, **kwargs_mapped)
            else:
                return self._astromodel_function.evaluate(
                    energy, **kwargs_mapped
                ).reshape(shape)

    @property
    def mapping(self):
        return self._mapping

    @property
    def mapping_free(self):
        return self._mapping_free
