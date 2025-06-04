import astropy.units as u
import numpy as np
from astromodels.core.units import get_units
from astromodels.sources import PointSource
from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,
    LogParabolaSpectralModel,
)
from scipy.integrate import dblquad

from gammapy_plugin.utils.astromodels_functions import (
    Exp_cutoff_powerlaw_gammapy,
    Gaussian_on_sphere,
    Log_parabola_gammapy,
)

get_units().energy = u.TeV


def test_log_parabola_gammapy():
    logp = Log_parabola_gammapy()
    # needed for using units
    ps = PointSource(source_name="test_source", spectral_shape=logp, ra=0, dec=0)
    logp.K = 1e-11 * u.Unit("TeV-1 cm-2 s-1")
    logp.piv = 2 * u.TeV
    logp.alpha.value = 2.0
    logp.beta.value = 0.9

    logp_ref = LogParabolaSpectralModel(
        amplitude=1e-11 * u.Unit("TeV-1 cm-2 s-1"),
        reference=2 * u.TeV,
        alpha=2.0,
        beta=0.9,
    )
    ref = np.geomspace(0.1, 100, 200) * u.TeV
    # assert this is also close when using no units in astromodels
    assert np.allclose(logp_ref(ref).value, logp(ref.value))
    assert np.isclose(ps.spectrum.main.Log_parabola_gammapy.K.value, 1e-11, atol=1e-20)


def test_exp_cutoff_powerlaw_gammapy():
    expc = Exp_cutoff_powerlaw_gammapy()
    # needed for using units
    ps = PointSource(source_name="test_source", spectral_shape=expc, ra=0, dec=0)
    expc.K = 1e-11 * u.Unit("TeV-1 cm-2 s-1")
    expc.piv = 2 * u.TeV
    expc.alpha.value = 0.89
    expc.lambda_ = 10 * 1 / u.TeV

    expc_ref = ExpCutoffPowerLawSpectralModel(
        amplitude=1e-11 * u.Unit("TeV-1 s-1 cm-2"),
        alpha=0.89,
        lambda_=10 * u.Unit("TeV-1"),
        reference=2 * u.TeV,
    )
    ref = np.geomspace(0.1, 100, 200) * u.TeV
    # assert this is also close when using no units in astromodels
    assert np.allclose(expc_ref(ref).value, expc(ref.value))
    assert np.isclose(
        ps.spectrum.main.Exp_cutoff_powerlaw_gammapy.K.value, 1e-11, atol=1e-20
    )


def test_gaussian_on_sphere():
    gauss = Gaussian_on_sphere(lon0=0, lat0=0, sigma=1)
    tot_int = gauss.get_total_spatial_integral(binsz=(1e-1, 1e-1))
    tot_int_ref, _ = dblquad(gauss, -20, 20, -20, 20, args={"sigma": 1})
    assert np.isclose(tot_int, tot_int_ref, rtol=1e-3)
