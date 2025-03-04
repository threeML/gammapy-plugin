import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion
from gammapy.data import DataStore
from gammapy.datasets import (
    Datasets,
    SpectrumDataset,
)
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import SkyModel
from gammapy_plugin.GammapyLike import GammapyLike
from gammapy.modeling.models import LogParabolaSpectralModel

# from threeML.bayesian.bayesian_analysis import BayesianAnalysis
from threeML.analysis_results import BayesianResults, MLEResults
from threeML.classicMLE.joint_likelihood import JointLikelihood
from threeML.data_list import DataList
from astromodels.core.model import Model
from astromodels.sources import PointSource
from astromodels.functions import Log_parabola, Log_uniform_prior, Uniform_prior


def get_close(threeML_results, gammapy_result_dict):
    if isinstance(threeML_results, BayesianResults):
        bm = threeML_results.get_median_fit_model().free_parameters
    elif isinstance(threeML_results, MLEResults):
        bm = threeML_results.optimized_model.free_parameters
    else:
        raise NotImplementedError
    hdp = threeML_results.get_data_frame(error_type="hpd")
    for p in bm.keys():
        pn = p.split(".")[-1]
        if pn == "K":
            pn = "amplitude"
        for gp in gammapy_result_dict["spectral"]["parameters"]:
            if gp["name"] == pn:
                break
        if pn == "amplitude":
            v = (bm[p].value * bm[p].unit).to("TeV-1 cm-2 s-1").value
            min_v = (
                ((hdp.loc[p]["negative_error"] + bm[p].value) * bm[p].unit)
                .to("TeV-1 cm-2 s-1")
                .value
            )
            max_v = (
                ((hdp.loc[p]["positive_error"] + bm[p].value) * bm[p].unit)
                .to("TeV-1 cm-2 s-1")
                .value
            )
        elif pn == "alpha":
            v = -bm[p].value
            min_v = hdp.loc[p]["negative_error"] - bm[p].value
            max_v = hdp.loc[p]["positive_error"] - bm[p].value
        else:
            v = bm[p].value
            min_v = hdp.loc[p]["negative_error"] + bm[p].value
            max_v = hdp.loc[p]["positive_error"] + bm[p].value
        return bool(gp["value"] <= max_v and min_v <= gp["value"])


def test_crab_spectrum():
    datastore = DataStore.from_dir("$GAMMAPY_DATA/hess-dl3-dr1/")
    obs_ids = [23523, 23526, 23559, 23592]
    observations = datastore.get_observations(obs_ids)
    target_position = SkyCoord(ra=83.63, dec=22.01, unit="deg", frame="icrs")

    on_region_radius = Angle("0.11 deg")
    on_region = CircleSkyRegion(center=target_position, radius=on_region_radius)
    exclusion_region = CircleSkyRegion(
        center=SkyCoord(183.604, -8.708, unit="deg", frame="galactic"),
        radius=0.5 * u.deg,
    )

    skydir = target_position.galactic
    geom = WcsGeom.create(
        npix=(150, 150), binsz=0.05, skydir=skydir, proj="TAN", frame="icrs"
    )

    exclusion_mask = ~geom.region_mask([exclusion_region])
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
    bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)
    safe_mask_maker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=15)
    datasets = Datasets()

    for obs_id, observation in zip(obs_ids, observations):
        dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
        dataset_on_off = bkg_maker.run(dataset, observation)
        dataset_on_off = safe_mask_maker.run(dataset_on_off, observation)
        datasets.append(dataset_on_off)

    logp = Log_parabola()
    logp.K.prior = Log_uniform_prior(lower_bound=1e-22, upper_bound=1e-19)
    logp.K.value = 1e-21
    logp.piv.value = 1e9
    logp.piv.free = False
    logp.alpha.prior = Uniform_prior(lower_bound=-3, upper_bound=-1)
    logp.alpha.value = -2
    logp.beta.prior = Uniform_prior(lower_bound=-2, upper_bound=2)
    logp.beta.value = 1
    ps = PointSource(
        source_name="crab",
        ra=target_position.ra.deg,
        dec=target_position.dec.deg,
        spectral_shape=logp,
    )
    model = Model(ps)
    gl = GammapyLike("hess", sources="crab")
    gl.set_datasets(datasets)
    gl.set_model(model)

    # ba = BayesianAnalysis(likelihood_model=model, data_list=DataList(gl))
    # ba.set_sampler("multinest")
    # ba.sampler.setup()
    # ba.sample()
    # res = ba.results

    jl = JointLikelihood(model, DataList(gl))
    jl.fit()
    res = jl.results

    logp_gammapy = LogParabolaSpectralModel(
        amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"), reference=1 * u.TeV
    )
    models = SkyModel(spectral_model=logp_gammapy, name="crab_gp")
    dataset_stacked = Datasets(datasets).stack_reduce()
    dataset_stacked.models = models
    fit_stacked = Fit()
    fit_stacked.run([dataset_stacked])

    return get_close(res, models.spectral_model.to_dict())
