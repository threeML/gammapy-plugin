import astropy.units as u
import numpy as np
import pytest
from astromodels.core.sky_direction import SkyDirection
from astromodels.functions import Gaussian_on_sphere, Powerlaw

from gammapy_plugin.models import (
    PointSourceModelConverted,
    SpatialModelConverted,
    SpectralModelConverted,
)


def test_spectral_model_converted():
    pl = Powerlaw()
    with pytest.raises(ValueError):
        SpectralModelConverted(pl)
    with pytest.raises(ValueError):
        spec = SpectralModelConverted(pl)
    pl.set_units(u.TeV, u.Unit("TeV-1 cm-2 s-1"))
    spec = SpectralModelConverted(pl)
    x = 0.5 * u.TeV
    para_dict = {}
    for v in pl.parameters.values():
        para_dict[v.path] = v
    assert spec.evaluate(x, **para_dict) == pl(x)
    x = np.ones(10).reshape((2, 5)) * 0.5 * u.TeV
    assert np.all(spec.evaluate(x, **para_dict) == pl(x))

    pl1 = Powerlaw()
    pl2 = Powerlaw()
    pl1.set_units(u.TeV, u.Unit("TeV-1 cm-2 s-1"))
    pl2.set_units(u.TeV, u.Unit("TeV-1 cm-2 s-1"))

    spec_comps = SpectralModelConverted([pl1, pl2])
    x = 0.5 * u.TeV
    para_dict = {}
    for pl_x in [pl1, pl2]:
        for v in pl_x.parameters.values():
            para_dict[v.path] = v
    assert spec_comps.evaluate(x, **para_dict) == pl1(x) + pl2(x)
    x = np.ones(10).reshape((2, 5)) * 0.5 * u.TeV
    assert np.all(spec_comps.evaluate(x, **para_dict) == pl1(x) + pl2(x))

    with pytest.raises(ValueError):
        pl2 = Powerlaw()
        spec_comps = SpectralModelConverted([pl1, pl2])

    with pytest.raises(AssertionError):
        pl2 = Powerlaw()
        pl2.set_units(u.TeV, None)
        spec_comps = SpectralModelConverted([pl1, pl2])
    with pytest.raises(NotImplementedError):
        string = "a string"
        SpectralModelConverted(string)


def test_spectral_model_integral_preserves_integral_units():
    pl = Powerlaw()
    pl.set_units(u.keV, u.Unit("keV-1 cm-2 s-1"))
    spec = SpectralModelConverted(pl)

    energy_min = np.array([1.0, 2.0]) * u.TeV
    energy_max = np.array([2.0, 4.0]) * u.TeV
    result = spec.integral(energy_min, energy_max)

    expected = u.Quantity(
        [
            pl.integrate(lo, hi).to_value(spec._integral_unit)
            for lo, hi in zip(energy_min, energy_max, strict=True)
        ],
        spec._integral_unit,
    )
    assert result.unit == spec._integral_unit
    assert np.allclose(result.value, expected.value)

    scalar_result = spec.integral(1 * u.TeV, 2 * u.TeV)
    assert scalar_result.isscalar
    assert scalar_result.unit == spec._integral_unit
    assert np.isclose(scalar_result.value, expected[0].value)


def test_spectral_model_integral_multiple_components():
    components = [Powerlaw(), Powerlaw()]
    for component in components:
        component.set_units(u.TeV, u.Unit("TeV-1 cm-2 s-1"))
    spec = SpectralModelConverted(components)

    actual = spec.integral(1 * u.TeV, 2 * u.TeV)
    expected = sum(
        (component.integrate(1 * u.TeV, 2 * u.TeV) for component in components),
        0 * spec._integral_unit,
    )

    assert actual.unit == spec._integral_unit
    assert np.isclose(actual.value, expected.to_value(spec._integral_unit))


def test_spatial_model_converted():
    gauss = Gaussian_on_sphere()
    from astromodels.core.parameter_transformation import LogarithmicTransformation

    gauss.sigma.value = 1.0
    gauss.sigma._transformation = LogarithmicTransformation()
    SpatialModelConverted(gauss)
    # FIXME:


def test_point_spatial_model_converted():
    sd = SkyDirection(ra=360, dec=0)
    sd.parameters["ra"].free = True
    PointSourceModelConverted(sd, frame="icrs")
    with pytest.raises(NotImplementedError):
        PointSourceModelConverted(sd, frame="fk5")
