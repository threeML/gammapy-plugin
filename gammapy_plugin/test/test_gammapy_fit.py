import astropy.units as u
import numpy as np
import pytest
from astromodels.core.model import Model
from astromodels.functions import Powerlaw
from astromodels.sources import PointSource
from gammapy.modeling import Fit
from gammapy.modeling.models import PowerLawSpectralModel, SkyModel

from gammapy_plugin.converter import AstromodelConverter


def test_gammapy_fit(crab_test_data):
    datasets = crab_test_data.copy()
    datasets_gp = crab_test_data.copy()
    np.random.seed(1234)

    pl = Powerlaw()
    ps = PointSource("crab", spectral_shape=pl, ra=83.63, dec=22.01)
    model = Model(ps)
    conv = AstromodelConverter(model)

    datasets.models = [conv.gammapy_models[0]]  # using ReflectedRegionsBackground

    fit = Fit()
    result = fit.run(datasets=datasets)

    spectral_model = PowerLawSpectralModel(
        amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
        index=2,
        reference=1 * u.TeV,
    )
    model_gp = SkyModel(spectral_model=spectral_model, name="crab")

    datasets_gp.models = [model_gp]

    fit_joint = Fit()
    result_joint = fit_joint.run(datasets=datasets_gp)

    assert np.isclose(
        result.models[0].parameters["crab.spectrum.main.Powerlaw.index"].value,
        -result_joint.models[0].parameters["index"].value,
        rtol=0.1,
    )


def test_xspec_wrapping(crab_test_data):

    pytest.importorskip("xspec")

    from astromodels.xspec import XS_powerlaw

    datasets = crab_test_data.copy()
    datasets = datasets.stack_reduce()
    datasets_gp = crab_test_data.copy()
    datasets_gp = datasets_gp.stack_reduce()
    np.random.seed(1234)
    pl = XS_powerlaw()
    ps = PointSource("crab", spectral_shape=pl, ra=83.633, dec=22.014)
    pl.phoindex = 2
    pl.phoindex.min_value = -np.nan
    pl.phoindex.max_value = np.nan

    model = Model(ps)
    conv = AstromodelConverter(model)
    datasets.models = [conv.gammapy_models[0]]  # using ReflectedRegionsBackground
    fit = Fit()
    result = fit.run(datasets=datasets)

    spectral_model = PowerLawSpectralModel(
        amplitude=100 * u.Unit("cm-2 s-1 keV-1"),
        index=2,
        reference=1 * u.keV,
    )
    model_gp = SkyModel(spectral_model=spectral_model, name="crab")

    datasets_gp.models = [model_gp]

    fit_joint = Fit()
    result_joint = fit_joint.run(datasets=datasets_gp)
    assert np.isclose(
        result.models[0].parameters["crab.spectrum.main.XS_powerlaw.norm"].value,
        result_joint.models[0].parameters["amplitude"].value,
        atol=result_joint.models[0].parameters["amplitude"].error,
    )
    assert np.isclose(
        result.models[0].parameters["crab.spectrum.main.XS_powerlaw.phoindex"].value,
        result_joint.models[0].parameters["index"].value,
        atol=result_joint.models[0].parameters["index"].error,
    )
