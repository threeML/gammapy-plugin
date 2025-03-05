import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion
from gammapy.data import DataStore
from gammapy.datasets import Datasets, SpectrumDataset
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import PowerLawSpectralModel, SkyModel
from gammapy_plugin.GammapyLike import GammapyLike
from astromodels.functions import (
    Powerlaw,
    Uniform_prior,
    Log_uniform_prior,
    Disk_on_sphere,
)
from astromodels.core.model import Model
from astromodels.sources.extended_source import ExtendedSource
from threeML import JointLikelihood, DataList
from gammapy_plugin.test.utils import get_close


def test_extended_source_spectrum():
    """
    Here we only fit the spectrum of the extended source RX J1713
    """
    datastore = DataStore.from_dir("$GAMMAPY_DATA/hess-dl3-dr1/")
    obs_ids = [20326, 20327, 20349, 20350, 20396, 20397]
    # In case you want to use all RX J1713 data in the H.E.S.S. DR1
    # other_ids=[20421, 20422, 20517, 20518, 20519, 20521, 20898, 20899, 20900]

    observations = datastore.get_observations(obs_ids)
    target_position = SkyCoord(347.3, -0.5, unit="deg", frame="galactic")
    radius = Angle("0.5 deg")
    on_region = CircleSkyRegion(target_position, radius)
    # The binning of the final spectrum is defined here.
    energy_axis = MapAxis.from_energy_bounds(0.1, 40.0, 10, unit="TeV")

    # Reduced IRFs are defined in true energy (i.e. not measured energy).
    energy_axis_true = MapAxis.from_energy_bounds(
        0.05, 100, 30, unit="TeV", name="energy_true"
    )

    geom = RegionGeom(on_region, axes=[energy_axis])

    dataset_empty = SpectrumDataset.create(
        geom=geom,
        energy_axis_true=energy_axis_true,
    )
    maker = SpectrumDatasetMaker(
        selection=["counts", "exposure", "edisp"], use_region_center=False
    )
    bkg_maker = ReflectedRegionsBackgroundMaker()
    safe_mask_maker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=10)
    datasets = Datasets()

    for obs in observations:
        # A SpectrumDataset is filled in this geometry
        dataset = maker.run(dataset_empty.copy(name=f"obs-{obs.obs_id}"), obs)

        # Define safe mask
        dataset = safe_mask_maker.run(dataset, obs)

        # Compute OFF
        dataset = bkg_maker.run(dataset, obs)

        # Append dataset to the list
        datasets.append(dataset)

    datasets_copy = datasets.copy()

    spectral_model = PowerLawSpectralModel(
        index=2, amplitude=2e-11 * u.Unit("cm-2 s-1 TeV-1"), reference=1 * u.TeV
    )
    model = SkyModel(spectral_model=spectral_model, name="RXJ 1713")

    datasets.models = [model]

    fit_joint = Fit()
    fit_joint.run(datasets=datasets)

    gl = GammapyLike("hess")
    gl.set_datasets(datasets_copy)
    pl = Powerlaw()
    pl.index.prior = Uniform_prior(lower_bound=-4, upper_bound=-1)
    pl.index.value = -3
    pl.K.prior = Log_uniform_prior(lower_bound=1e-22, upper_bound=1e-19)
    pl.K.value = 1e-21
    pl.piv.value = 1e9
    pl.piv.free = False
    disk = Disk_on_sphere(
        lon0=target_position.transform_to("icrs").ra.deg,
        lat0=target_position.transform_to("icrs").dec.deg,
        radius=radius.value,
    )
    disk.lon0.free = False
    disk.lat0.free = False
    disk.radius.free = False
    es = ExtendedSource(source_name="rxj1713", spatial_shape=disk, spectral_shape=pl)
    model_am = Model(es)
    gl.set_model(model_am)
    jl = JointLikelihood(model_am, DataList(gl))
    jl.fit()
    res = jl.results

    assert get_close(res, model.spectral_model.to_dict()) is True
