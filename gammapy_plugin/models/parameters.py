import numpy as np
from gammapy.modeling.parameter import Parameter, Parameters


class ParameterWrapped(Parameter):
    def __init__(self, *args, **kwargs):
        self._am_model = kwargs.pop("astromodels_model", None)
        self._am_factor = kwargs.pop("astromodels_factor", 1.0)
        super().__init__(*args, **kwargs)

    @property
    def factor(self):
        """Factor as a float."""
        return self._factor

    @factor.setter
    def factor(self, val):
        self._factor = float(np.nan_to_num(val))
        self._value = float(self.inverse_transform(self._factor))
        if self._am_model is not None:
            # print(
            #   f"Setting new value via factor {self._factor}: {self._value/self._am_factor} for {self._name}"
            # )

            self._am_model.parameters[self.name] = self._value / self._am_factor

    @property
    def value(self):
        """Value = factor x scale (float)."""
        return self._value

    @value.setter
    def value(self, val):
        self._value = float(val)
        self._factor = self.transform(val)
        if self._am_model is not None:
            # print(f"Setting new value {self._value} for {self._name}")
            self._am_model.parameters[self.name] = self._value / self._am_factor


class ParametersWrapped(Parameters):
    """
    Class for updating the astromodels parameters according to the gammapy parameters.
    """

    def set_parameter_factors(self, factors):
        """Set factor of all parameters.

        Used in the optimizer interface.
        """
        idx = 0
        for parameter in self._parameters:
            if not parameter.frozen:
                parameter.factor = factors[idx]
                idx += 1
