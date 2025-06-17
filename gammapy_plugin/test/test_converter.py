from astromodels.core.model import Model
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
