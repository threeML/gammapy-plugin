import pytest
from astromodels.core.model import Model
from astromodels.functions import Powerlaw
from astromodels.sources import PointSource
from gammapy.modeling.models import FoVBackgroundModel, SkyModel

from gammapy_plugin.converter import AstromodelConverter
from gammapy_plugin.GammapyLike import GammapyLike
from gammapy_plugin.test.utils import read_in_gammapy_datasets
from gammapy_plugin.utils.package_data import get_path_of_data_dir


def test_set_datasets_stacked():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("test/rxj17137_3946/")
    )
    gl_stacked = GammapyLike(name="stacked")
    gl_stacked.set_datasets(datasets.stack_reduce(name="stacked"), mode="stacked")
    gl_stacked.set_datasets(datasets.stack_reduce(name="stacked"), mode="individual")
    gl_stacked.set_datasets(datasets, mode="stacked")


def test_set_datasets_individual():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("test/rxj17137_3946/")
    )
    gl_individual = GammapyLike(name="individual")
    gl_individual.set_datasets(datasets)


def test_set_datasets_list():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("test/rxj17137_3946/")
    )
    gl_list = GammapyLike(name="list")
    datasets_list = [d for d in datasets]
    gl_list.set_datasets(datasets_list, mode="stacked")


def test_set_multiple_datasets():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("test/rxj17137_3946/")
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
    gl_multi_individual._update_background_models()


def test_set_datasets_error():
    with pytest.raises(TypeError):
        gl_error = GammapyLike(name="error")
        gl_error.set_datasets(None)


def test_set_datasets_wrong_mode():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("test/rxj17137_3946/")
    )
    gl_wrong_mode = GammapyLike("test")
    with pytest.raises(ValueError):
        gl_wrong_mode.set_datasets(datasets, mode="something_weird")


def test_set_sources():
    gl = GammapyLike(name="test")
    gl.set_sources()
    gl.set_sources(["test"])
    gl.set_sources("test")
    with pytest.raises(ValueError):
        gl.set_sources(GammapyLike("hehe"))


def test_set_background_models_init():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("test/rxj17137_3946/")
    )
    bkg_models = {}
    for d in datasets:
        bkg_models[d.name + "_bkg"] = FoVBackgroundModel(
            name=d.name + "_bkg", dataset_name=d.name
        )
    GammapyLike(name="init_bkg", background_models=bkg_models)


def test_wrong_bkg_model():
    gl = GammapyLike(name="wrong_bkg")
    with pytest.raises(TypeError):
        gl.set_background_models("string")


def test_properties():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("test/rxj17137_3946/")
    )
    bkg_models = {}
    for d in datasets:
        bkg_models[d.name + "_bkg"] = FoVBackgroundModel(
            name=d.name + "_bkg", dataset_name=d.name
        )
    gl = GammapyLike(name="init_bkg")
    gl.set_datasets(datasets, mode="individual")
    gl.gammapy_model
    gl.set_background_models(list(bkg_models.values()))
    pl = Powerlaw()
    ps = PointSource(source_name="test", ra=0, dec=0, spectral_shape=pl)
    model = Model(ps)
    conv = AstromodelConverter(model)
    gl.set_model(model, conv)
    assert gl.model == model
    assert gl.astromodel_converter == conv
    assert gl.frame == "icrs"
    assert isinstance(gl.gammapy_model[0], SkyModel)
