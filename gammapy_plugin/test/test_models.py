import astropy.units as u
import pytest
from astromodels.core.sky_direction import SkyDirection
from astromodels.functions import Powerlaw, Gaussian_on_sphere

from gammapy_plugin.models import (
    PointSourceModelConverted,
    SpatialModelConverted,
    SpectralModelConverted,
    TemporalModelConverted,
)


def test_spectral_model_converted():
    pl = Powerlaw()
    with pytest.raises(ValueError):
        SpectralModelConverted(pl)
    pl.set_units(u.TeV, u.Unit("TeV-1 cm-2 s-1"))
    spec = SpectralModelConverted(pl)
    x = 0.5 * u.TeV
    assert spec.evaluate(x, **{"K": pl.K, "index": pl.index, "piv": pl.piv}) == pl(x)


def test_spatial_model_converted():
    gauss = Gaussian_on_sphere()
    from astromodels.core.parameter_transformation import LogarithmicTransformation

    gauss.sigma._transformation = LogarithmicTransformation()
    SpatialModelConverted(gauss)


def test_point_spatial_model_converted():
    sd = SkyDirection(ra=360, dec=0)
    sd.parameters["ra"].free = True
    PointSourceModelConverted(sd, frame="icrs")
    with pytest.raises(NotImplementedError):
        PointSourceModelConverted(sd, frame="fk5")


def test_temporal_model_converted():
    with pytest.raises(NotImplementedError):
        pl = Powerlaw()
        TemporalModelConverted(pl)
