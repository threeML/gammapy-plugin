import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from gammapy.data.data_store import DataStore
from gammapy.datasets import Datasets, MapDataset
from gammapy.makers import FoVBackgroundMaker, MapDatasetMaker, SafeMaskMaker
from gammapy.maps import MapAxis
from gammapy.maps.wcs.geom import WcsGeom
from gammapy.modeling.models import FoVBackgroundModel
from regions import CircleSkyRegion

from gammapy_plugin.GammapyLike import GammapyLike
from gammapy_plugin.test.utils import read_in_gammapy_datasets
from gammapy_plugin.utils.package_data import get_path_of_data_file

datasets = read_in_gammapy_datasets(
    get_path_of_data_file("datasets/test/rxj17137_3946/")
)


def test_set_datasets_stacked():
    gl_stacked = GammapyLike(name="stacked")
    gl_stacked.set_datasets(datasets.stack_reduce(name="stacked"), mode="stacked")
    gl_stacked.set_datasets(datasets.stack_reduce(name="stacked"), mode="individual")
    gl_stacked.set_datasets(datasets, mode="stacked")


def test_set_datasets_individual():
    gl_individual = GammapyLike(name="individual")
    gl_individual.set_datasets(datasets)
    gl_individual.set_datasets(datasets[0])


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
