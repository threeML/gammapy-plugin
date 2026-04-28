import astropy.units as u
import numpy as np
import pytest
from astromodels.core.model import Model
from astromodels.core.units import get_units
from astromodels.functions import Log_uniform_prior, Powerlaw, Uniform_prior
from astromodels.sources import PointSource
from astropy.coordinates import Angle, SkyCoord
from gammapy.data import DataStore
from gammapy.datasets import Datasets, SpectrumDataset
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from regions import CircleSkyRegion
from threeML import BayesianAnalysis
from threeML.data_list import DataList

from gammapy_plugin.GammapyLike import GammapyLike


@pytest.fixture
def _crab_dataset():
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


@pytest.fixture
def _pks2155304_dataset():
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


def test_crab_spectrum(_crab_dataset, tmp_path):
    get_units().energy = u.keV

    position = SkyCoord.from_name("Crab")
    pl = Powerlaw()
    ps = PointSource(
        "crab", ra=position.ra.deg, dec=position.dec.deg, spectral_shape=pl
    )
    pl.K.prior = Log_uniform_prior(lower_bound=1e-23, upper_bound=1e-19)
    pl.index.prior = Uniform_prior(lower_bound=-4, upper_bound=-0.05)
    pl.piv.value = 1.45 * 1e9
    pl.piv.free = False
    model = Model(ps)

    gl = GammapyLike(name="gammapy_plugin")
    gl.set_datasets(_crab_dataset, mode="stacked")
    gl.set_sources("crab")
    gl.set_model(model)

    dl = DataList(gl)
    ba = BayesianAnalysis(model, dl)
    ba.set_sampler("multinest")
    ba.sampler.setup(resume=False, chain_name=str(tmp_path / "chain/fit-"))
    ba.sample(quiet=True)

    assert np.isclose(
        ba.results.optimized_model.crab.spectrum.main.Powerlaw.index.value,
        -2.63,
        atol=0.07,
    )
    assert np.isclose(
        ba.results.optimized_model.crab.spectrum.main.Powerlaw.K.value,
        16.3 * 1e-21,
        atol=0.9 * 1e-21,
    )


def test_pks2155304_spectrum(_pks2155304_dataset, tmp_path):
    get_units().energy = u.keV

    position = SkyCoord.from_name("PKS 2155-304")
    pl = Powerlaw()
    ps = PointSource("pks", ra=position.ra.deg, dec=position.dec.deg, spectral_shape=pl)
    pl.K.prior = Log_uniform_prior(lower_bound=1e-25, upper_bound=1e-19)
    pl.index.prior = Uniform_prior(lower_bound=-7, upper_bound=-1.0)
    pl.piv.value = 0.65 * 1e9
    pl.piv.free = False
    model = Model(ps)
    for d in _pks2155304_dataset:
        print(d.energy_range)
    gl = GammapyLike(name="gammapy_plugin")
    gl.set_datasets(_pks2155304_dataset, mode="stacked")
    gl.set_sources("pks")
    gl.set_model(model)

    dl = DataList(gl)
    ba = BayesianAnalysis(model, dl)
    ba.set_sampler("multinest")
    ba.sampler.setup(resume=False, chain_name=str(tmp_path / "pks/fit-"))
    ba.sample(quiet=True)

    assert np.isclose(
        ba.results.optimized_model.pks.spectrum.main.Powerlaw.index.value,
        -3.46,
        atol=0.25,
    )
    assert np.isclose(
        ba.results.optimized_model.pks.spectrum.main.Powerlaw.K.value,
        17.5 * 1e-21,
        atol=2.1 * 1e-21,
    )
