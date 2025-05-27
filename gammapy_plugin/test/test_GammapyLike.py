from gammapy_plugin.GammapyLike import GammapyLike
from gammapy.data.data_store import DataStore
from gammapy.maps.wcs.geom import WcsGeom
from gammapy.maps import MapAxis
from astropy.coordinates import SkyCoord
import astropy.units as u
from gammapy.datasets import Datasets, MapDataset
from gammapy.modeling.models import FoVBackgroundModel
from gammapy.makers import FoVBackgroundMaker, MapDatasetMaker, SafeMaskMaker
from regions import CircleSkyRegion
import pytest

datastore = DataStore.from_dir("$GAMMAPY_DATA/hess-dl3-dr1/")
target_position = SkyCoord.from_name("RX J1713.7-3946").galactic

selection = dict(
    type="sky_circle",
    frame="galactic",
    lon=target_position.l,
    lat=target_position.b,
    radius="5deg",
)
select_obs_tab = datastore.obs_table.select_observations(selection)

obs = datastore.get_observations(select_obs_tab["OBS_ID"])

# Prepare the geometry
energy_axis = MapAxis.from_energy_bounds(0.3, 10.0, 15, unit="TeV")
energy_axis_true = MapAxis.from_energy_bounds(
    0.1, 20, 10, per_decade=True, unit="TeV", name="energy_true"
)
geom = WcsGeom.create(
    skydir=target_position,
    binsz=0.02,
    width=(6 * u.deg, 6 * u.deg),
    frame="galactic",
    axes=[energy_axis],
)
circle = CircleSkyRegion(center=target_position, radius=1 * u.deg)
regions = [circle]
exclusion_mask = ~geom.region_mask(regions=regions)
maker = MapDatasetMaker(
    selection=["counts", "background", "psf", "edisp", "exposure"],
)
safe_mask_maker = SafeMaskMaker(
    methods=["offset-max", "aeff-max", "bkg-peak"], offset_max="2.3 deg"
)
fov_bkg_maker = FoVBackgroundMaker(method="fit", exclusion_mask=exclusion_mask)

datasets = Datasets()
for o in obs:
    dataset = MapDataset.create(
        geom=geom, energy_axis_true=energy_axis_true, name=f"HESS_{o.obs_id}"
    )
    dataset = maker.run(dataset, o)
    dataset = safe_mask_maker.run(dataset, o)
    bkg_model = FoVBackgroundModel(name=f"{o.obs_id}_bkg", dataset_name=dataset.name)
    dataset.models = [bkg_model]
    dataset = fov_bkg_maker.run(dataset)
    datasets.append(dataset)


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
