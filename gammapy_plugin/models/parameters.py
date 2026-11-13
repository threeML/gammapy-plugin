import astropy.units as u
from gammapy.modeling.parameter import Parameter


class ParameterWrapped(Parameter):
    """Gammapy parameter linked to an astromodels parameter."""

    def __init__(self, *args, **kwargs):
        # Do not bind during construction:
        # Gammapy will deepcopy these parameters.
        self._am_parameter = None
        super().__init__(*args, **kwargs)

    def bind(self, am_parameter):
        """Connect to the actual astromodels Parameter object."""
        if not self.unit.is_equivalent(am_parameter.unit):
            raise u.UnitConversionError(
                f"Incompatible units: {self.unit} and {am_parameter.unit}"
            )

        self._am_parameter = am_parameter

    def _sync_to_astromodels(self):
        if self._am_parameter is None:
            return

        # Convert into the units used by astromodels.
        value = self.quantity.to_value(self._am_parameter.unit)
        self._am_parameter.value = value

    @Parameter.value.setter
    def value(self, val):
        Parameter.value.fset(self, val)
        self._sync_to_astromodels()

    @Parameter.factor.setter
    def factor(self, val):
        Parameter.factor.fset(self, val)
        self._sync_to_astromodels()

    @Parameter.quantity.setter
    def quantity(self, val):
        quantity = u.Quantity(val)

        if not quantity.unit.is_equivalent(self.unit):
            raise u.UnitConversionError(f"Incompatible unit for {self.name}")

        # Set unit before value, so synchronization
        # always uses the correct physical unit.
        self.unit = quantity.unit
        self.value = quantity.value
