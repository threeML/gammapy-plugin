import astropy.units as u
import numpy as np
from astromodels.core.model import Model
from astromodels.core.spectral_component import SpectralComponent
from astromodels.functions import Powerlaw
from astromodels.sources import PointSource

from gammapy_plugin.converter import AstromodelConverter


def test_source_converter():
    pl = Powerlaw()
    ps = PointSource(ra=0, dec=0, spectral_shape=pl, source_name="test_ps")
    model = Model(ps)
    conv = AstromodelConverter(model, frame="galactic")
    sc = conv._converted_sources["test_ps"]
    assert sc.astromodels_source == ps
    assert len(conv.converted_sources) > 0


def test_parameter_mapping():
    pl1 = Powerlaw(index=-9.9, K=12)
    pl2 = Powerlaw()
    pl2.K
    sc1 = SpectralComponent("comp1", pl1)
    sc2 = SpectralComponent("comp2", pl2)
    ps1 = PointSource(ra=0, dec=0, components=[sc1, sc2], source_name="test_ps")
    ps2 = PointSource(
        ra=10,
        dec=10,
        spectral_shape=Powerlaw(index=0, K=1),
        source_name="source2",
    )
    model = Model(*[ps1, ps2])
    conv = AstromodelConverter(model, frame="galactic")
    pl2.K.value = 100
    conv._update_parameters()
    assert (
        conv._converted_sources["test_ps"]
        .astromodels_source.components["comp2"]
        .shape.K.value
        == 100
    )
    pl1.index.value = 3.1415
    assert (
        conv._converted_sources["test_ps"]
        .astromodels_source.components["comp1"]
        .shape.index.value
        == 3.1415
    )
    assert (
        conv._converted_sources[
            "source2"
        ].astromodels_source.spectrum.main.shape.index.value
        == 0
    )


def test_multi_comp():
    pl1 = Powerlaw(index=-1, K=1e-9, piv=1e9)
    pl2 = Powerlaw(index=-2, K=2e2 - 9, piv=1e9)
    sc1 = SpectralComponent("comp1", pl1)
    sc2 = SpectralComponent("comp2", pl2)
    ps1 = PointSource(ra=0, dec=0, components=[sc1, sc2], source_name="test_ps")
    model = Model(ps1)
    conv = AstromodelConverter(model, frame="galactic")
    assert np.isclose(
        ps1(1e9),
        conv.gammapy_models[0].spectral_model(1 * u.TeV).value,
    )
