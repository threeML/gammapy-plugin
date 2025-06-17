import pytest
from astromodels.core.model import Model
from astromodels.functions import Powerlaw
from astromodels.sources import PointSource
from gammapy.modeling.models import FoVBackgroundModel

from gammapy_plugin.GammapyLike import GammapyLike
from gammapy_plugin.test.utils import read_in_gammapy_datasets
from gammapy_plugin.utils.package_data import get_path_of_data_dir


def test_set_datasets_stacked():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("datasets/test/rxj17137_3946/")
    )
    gl_stacked = GammapyLike(name="stacked")
    gl_stacked.set_datasets(datasets.stack_reduce(name="stacked"), mode="stacked")
    gl_stacked.set_datasets(datasets.stack_reduce(name="stacked"), mode="individual")
    gl_stacked.set_datasets(datasets, mode="stacked")


def test_set_datasets_individual():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("datasets/test/rxj17137_3946/")
    )
    gl_individual = GammapyLike(name="individual")
    gl_individual.set_datasets(datasets)


def test_set_multiple_datasets():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("datasets/test/rxj17137_3946/")
    )
    gl_multi_individual = GammapyLike(name="individual")
    gl_multi_individual.set_datasets(datasets)
    pl = Powerlaw()
    ps = PointSource(ra=0, dec=0, spectral_shape=pl, source_name="test")
    model = Model(ps)
    bkg_models = []
    for d in datasets:
        bkg_models.append(FoVBackgroundModel(name=f"{d.name}_bkg", dataset_name=d.name))
    gl_multi_individual.set_model(model)
    gl_multi_individual.set_background_models(bkg_models)
    for d in gl_multi_individual.datasets:
        assert len(d.models) == 2, f"Too many models in dataset {d.name}"


def test_set_datasets_error():
    with pytest.raises(TypeError):
        gl_error = GammapyLike(name="error")
        gl_error.set_datasets(None)


def test_set_sources():
    gl = GammapyLike(name="test")
    gl.set_sources()
    gl.set_sources(["test"])
    gl.set_sources("test")
    with pytest.raises(ValueError):
        gl.set_sources(GammapyLike("hehe"))
