import astropy.units as u
import numpy as np
from astromodels.core.model import Model
from astromodels.core.units import get_units
from astromodels.functions import Log_uniform_prior, Powerlaw, Uniform_prior
from astromodels.sources import PointSource
from astropy.coordinates import SkyCoord
from threeML import BayesianAnalysis
from threeML.data_list import DataList

from gammapy_plugin.GammapyLike import GammapyLike


def test_crab_spectrum(crab_test_data, tmp_path):
    get_units().energy = u.keV

    position = SkyCoord.from_name("Crab")
    pl = Powerlaw()
    ps = PointSource(
        "crab", ra=position.ra.deg, dec=position.dec.deg, spectral_shape=pl
    )
    pl.K.prior = Log_uniform_prior(lower_bound=1e-23, upper_bound=1e-19)
    pl.index.prior = Uniform_prior(lower_bound=-4, upper_bound=-0.05)
    pl.piv.value = 1.45 * 1e9
    pl.piv.free = False
    model = Model(ps)

    gl = GammapyLike(name="gammapy_plugin")
    gl.set_datasets(crab_test_data, mode="stacked")
    gl.set_sources("crab")
    gl.set_model(model)

    dl = DataList(gl)
    ba = BayesianAnalysis(model, dl)
    ba.set_sampler("multinest")
    ba.sampler.setup(resume=False, chain_name=str(tmp_path / "chain/fit-"))
    ba.sample(quiet=True)

    assert np.isclose(
        ba.results.optimized_model.crab.spectrum.main.Powerlaw.index.value,
        -2.63,
        atol=0.07,
    )
    assert np.isclose(
        ba.results.optimized_model.crab.spectrum.main.Powerlaw.K.value,
        16.3 * 1e-21,
        atol=0.9 * 1e-21,
    )


def test_pks2155304_spectrum(pks2155304_test_data, tmp_path):
    get_units().energy = u.keV

    position = SkyCoord.from_name("PKS 2155-304")
    pl = Powerlaw()
    ps = PointSource("pks", ra=position.ra.deg, dec=position.dec.deg, spectral_shape=pl)
    pl.K.prior = Log_uniform_prior(lower_bound=1e-25, upper_bound=1e-19)
    pl.index.prior = Uniform_prior(lower_bound=-7, upper_bound=-1.0)
    pl.piv.value = 0.65 * 1e9
    pl.piv.free = False
    model = Model(ps)
    gl = GammapyLike(name="gammapy_plugin")
    gl.set_datasets(pks2155304_test_data, mode="stacked")
    gl.set_sources("pks")
    gl.set_model(model)

    dl = DataList(gl)
    ba = BayesianAnalysis(model, dl)
    ba.set_sampler("multinest")
    ba.sampler.setup(resume=False, chain_name=str(tmp_path / "pks/fit-"))
    ba.sample(quiet=True)

    assert np.isclose(
        ba.results.optimized_model.pks.spectrum.main.Powerlaw.index.value,
        -3.46,
        atol=0.25,
    )
    assert np.isclose(
        ba.results.optimized_model.pks.spectrum.main.Powerlaw.K.value,
        17.5 * 1e-21,
        atol=2.1 * 1e-21,
    )
