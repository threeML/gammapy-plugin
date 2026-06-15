import logging

import astropy.units as u
import pytest
from astropy.coordinates import Angle, SkyCoord
from filelock import FileLock
from gammapy.data import DataStore
from gammapy.datasets import Datasets, MapDataset, SpectrumDataset
from gammapy.makers import (
    MapDatasetMaker,
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from regions import CircleSkyRegion
from filelock import FileLock

log = logging.getLogger(__name__)


def produce_rxj_test_data():
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
    maker = MapDatasetMaker(
        selection=["counts", "background", "psf", "edisp", "exposure"],
    )
    safe_mask_maker = SafeMaskMaker(
        methods=["offset-max", "aeff-max", "bkg-peak"], offset_max="2.3 deg"
    )

    datasets = Datasets()
    for o in obs:
        dataset = MapDataset.create(
            geom=geom, energy_axis_true=energy_axis_true, name=f"HESS_{o.obs_id}"
        )
        dataset = maker.run(dataset, o)
        dataset = safe_mask_maker.run(dataset, o)
        datasets.append(dataset)

    return datasets


@pytest.fixture(scope="session")
def rxj_test_data(tmp_path_factory, worker_id):
    if worker_id == "master":
        return produce_rxj_test_data()

    root_tmp_dir = tmp_path_factory.getbasetemp().parent
    fn = root_tmp_dir / "rxj.fits"
    with FileLock(str(fn) + ".lock"):
        if fn.is_file():
            data = Datasets.read(fn)
        else:
            data = produce_rxj_test_data()
            data.write(fn)
    return data


def produce_crab_test_data():
    datastore = DataStore.from_dir("$GAMMAPY_DATA/hess-dl3-dr1/")
    selection = dict(
        type="sky_circle",
        frame="icrs",
        lon="83.633 deg",
        lat="22.014 deg",
        radius="5 deg",
    )
    selected_obs_table = datastore.obs_table.select_observations(selection)
    observations = datastore.get_observations(selected_obs_table["OBS_ID"])
    target_position = SkyCoord(ra=83.63, dec=22.01, unit="deg", frame="icrs")

    on_region_radius = Angle("0.125 deg")
    on_region = CircleSkyRegion(
        center=target_position.galactic, radius=on_region_radius
    )

    skydir = target_position.galactic
    geom = WcsGeom.create(
        npix=(150, 150), binsz=0.05, skydir=skydir, proj="TAN", frame="galactic"
    )

    energy_axis = MapAxis.from_energy_bounds(
        0.5, 40, nbin=10, per_decade=True, unit="TeV", name="energy"
    )
    energy_axis_true = MapAxis.from_energy_bounds(
        0.1, 100, nbin=20, per_decade=True, unit="TeV", name="energy_true"
    )

    geom = RegionGeom.create(region=on_region, axes=[energy_axis])
    dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

    dataset_maker = SpectrumDatasetMaker(
        containment_correction=True, selection=["counts", "exposure", "edisp"]
    )
    bkg_maker = ReflectedRegionsBackgroundMaker()
    safe_mask_maker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=15)
    datasets = Datasets()

    for obs_id, observation in zip(selected_obs_table["OBS_ID"], observations):
        dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
        dataset_on_off = bkg_maker.run(dataset, observation)
        dataset_on_off = safe_mask_maker.run(dataset_on_off, observation)
        datasets.append(dataset_on_off)

    return datasets


@pytest.fixture(scope="session")
def crab_test_data(tmp_path_factory, worker_id):
    if worker_id == "master":
        return produce_crab_test_data()

    root_tmp_dir = tmp_path_factory.getbasetemp().parent
    fn = root_tmp_dir / "crab.fits"
    with FileLock(str(fn) + ".lock"):
        if fn.is_file():
            data = Datasets.read(fn)
        else:
            data = produce_crab_test_data()
            data.write(fn)
    return data


def produce_pks2155304_test_data():
    datastore = DataStore.from_dir("$GAMMAPY_DATA/hess-dl3-dr1/")
    selection = dict(
        type="sky_circle",
        frame="icrs",
        lon="329.72 deg",
        lat="-30.22 deg",
        radius="5 deg",
    )
    selected_obs_table = datastore.obs_table.select_observations(selection)
    selected_obs_table = selected_obs_table[
        selected_obs_table["TARGET_NAME"] == "PKS 2155-304 (steady)"
    ]

    observations = datastore.get_observations(selected_obs_table["OBS_ID"])
    target_position = SkyCoord.from_name("PKS 2155-304")

    on_region_radius = Angle("0.125 deg")
    on_region = CircleSkyRegion(
        center=target_position.galactic, radius=on_region_radius
    )

    skydir = target_position.galactic
    geom = WcsGeom.create(
        npix=(150, 150), binsz=0.05, skydir=skydir, proj="TAN", frame="galactic"
    )

    energy_axis = MapAxis.from_energy_bounds(
        0.4, 7.5, nbin=10, unit="TeV", name="energy"
    )
    energy_axis_true = MapAxis.from_energy_bounds(
        0.1, 100, nbin=20, per_decade=True, unit="TeV", name="energy_true"
    )

    geom = RegionGeom.create(region=on_region, axes=[energy_axis])
    dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

    dataset_maker = SpectrumDatasetMaker(
        containment_correction=True, selection=["counts", "exposure", "edisp"]
    )
    bkg_maker = ReflectedRegionsBackgroundMaker()
    safe_mask_maker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=15)
    datasets = Datasets()

    for obs_id, observation in zip(selected_obs_table["OBS_ID"], observations):
        dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
        dataset_on_off = bkg_maker.run(dataset, observation)
        dataset_on_off = safe_mask_maker.run(dataset_on_off, observation)
        datasets.append(dataset_on_off)

    return datasets


@pytest.fixture(scope="session")
def pks2155304_test_data(tmp_path_factory, worker_id):
    if worker_id == "master":
        return produce_pks2155304_test_data()

    root_tmp_dir = tmp_path_factory.getbasetemp().parent
    fn = root_tmp_dir / "pks2155304.fits"
    with FileLock(str(fn) + ".lock"):
        if fn.is_file():
            data = Datasets.read(fn)
        else:
            data = produce_pks2155304_test_data()
            data.write(fn)
    return data
