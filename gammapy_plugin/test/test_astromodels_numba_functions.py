import astropy.units as u
import numpy as np
from astromodels.core.model import Model
from astromodels.functions import Cutoff_powerlaw
from astromodels.sources import PointSource

from gammapy_plugin.converter import AstromodelConverter
from gammapy_plugin.GammapyLike import GammapyLike
from gammapy_plugin.test.utils import read_in_gammapy_datasets
from gammapy_plugin.utils.package_data import get_path_of_data_dir


def test_astromodels_numba_function():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("datasets/test/rxj17137_3946/")
    )
    gl = GammapyLike("test")
    gl.set_datasets(datasets[0])
    cpl = Cutoff_powerlaw()
    ps = PointSource(source_name="test", ra=0, dec=0, spectral_shape=cpl)
    model = Model(ps)
    conv = AstromodelConverter(model, convert_ps=False)
    gl.set_model(model, conv)
    res = (
        gl.datasets[0]
        .models[0]
        .evaluate(
            0 * u.deg, 0 * u.deg, np.geomspace(0.1, 40, 100).reshape(20, 5) * u.TeV
        )
    )
    assert res.shape == (20, 5)
