import os

import astropy.units as u
import numpy as np
from astromodels.core.model import Model
from astromodels.core.units import get_units
from astromodels.functions import Log_uniform_prior, Powerlaw, Uniform_prior
from astromodels.sources.extended_source import ExtendedSource
from astropy.coordinates import SkyCoord
from gammapy.makers import FoVBackgroundMaker
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    FoVBackgroundModel,
    GaussianSpatialModel,
    PowerLawSpectralModel,
    SkyModel,
)
from regions import CircleSkyRegion
from threeML import DataList, JointLikelihood

from gammapy_plugin.converter import AstromodelConverter
from gammapy_plugin.GammapyLike import GammapyLike
from gammapy_plugin.test.utils import read_in_gammapy_datasets
from gammapy_plugin.utils.astromodels_functions import Gaussian_on_sphere
from gammapy_plugin.utils.package_data import get_path_of_data_dir

get_units().energy = u.TeV
target_position = SkyCoord.from_name("RX J1713.7-3946").galactic


def test_extended_source_no_fov_bkg():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("datasets/test/rxj17137_3946/")
    )

    geom = datasets[0].geoms["geom"]
    circle = CircleSkyRegion(center=target_position, radius=1 * u.deg)
    regions = [circle]
    exclusion_mask = ~geom.region_mask(regions=regions)
    fov_bkg_maker = FoVBackgroundMaker(method="fit", exclusion_mask=exclusion_mask)
    for dataset in datasets:
        dataset = fov_bkg_maker.run(dataset)
    stacked = datasets.stack_reduce(name="stacked")
    pl = Powerlaw()
    spat = Gaussian_on_sphere(
        lon0=target_position.transform_to("galactic").l.deg,
        lat0=target_position.transform_to("galactic").b.deg,
        sigma=0.25,
    )
    es = ExtendedSource(source_name="rxj1713", spectral_shape=pl, spatial_shape=spat)
    pl.index.value = -2
    pl.index.prior = Uniform_prior(lower_bound=-3, upper_bound=-1)
    pl.K = 2.3 * 1e-11
    pl.K.prior = Log_uniform_prior(lower_bound=1e-18, upper_bound=1e-8)
    pl.piv.value = 1
    pl.piv.free = False
    spat.lon0.free = False
    spat.lat0.free = False
    spat.sigma.free = True
    spat.sigma.prior = Log_uniform_prior(lower_bound=0.1, upper_bound=1.0)
    model = Model(es)
    gl = GammapyLike("hess", frame="galactic")
    gl.set_datasets(
        stacked.copy(), mode="stacked"
    )  # making a copy to not interfer with the gammpy fit later
    gl.set_sources("rxj1713")
    conv = AstromodelConverter(model, frame="galactic")
    gl.set_model(model, conv)
    jl = JointLikelihood(model, DataList(gl))
    jl.fit()
    res = jl.results

    pl_gp = PowerLawSpectralModel(reference=1 * u.TeV)
    gauss_gp = GaussianSpatialModel(
        lon_0=347.269 * u.deg,
        lat_0=-0.257 * u.deg,
        frame="galactic",
    )
    gauss_gp.e.frozen = True
    gauss_gp.phi.frozen = True
    gauss_gp.lon_0.frozen = True
    gauss_gp.lat_0.frozen = True

    model_gp = SkyModel(name="rxj_gp", spatial_model=gauss_gp, spectral_model=pl_gp)
    ds_gp = stacked.copy()
    ds_gp.models = [model_gp]
    fit = Fit()
    resu = fit.run(datasets=ds_gp)

    assert (
        np.isclose(
            res.optimized_model.free_parameters[
                "rxj1713.Gaussian_on_sphere.sigma"
            ].value,
            resu.models.parameters["sigma"].value,
            rtol=1e-3,
        )
        is np.True_
    )
    assert (
        np.isclose(
            -res.optimized_model.free_parameters[
                "rxj1713.spectrum.main.Powerlaw.index"
            ].value,
            resu.models.parameters["index"].value,
            rtol=1e-3,
        )
        is np.True_
    )
    assert (
        np.isclose(
            res.optimized_model.free_parameters[
                "rxj1713.spectrum.main.Powerlaw.K"
            ].value
            * res.optimized_model.extended_sources[
                "rxj1713"
            ].spatial_shape.get_total_spatial_integral(1),
            resu.models.parameters["amplitude"].value,
            rtol=1e-3,
            atol=1e-20,
        )
        is np.True_
    )


def test_fov_bkg_model_setting():
    datasets = read_in_gammapy_datasets(
        get_path_of_data_dir().joinpath("datasets/test/rxj17137_3946/")
    )
    geom = datasets[0].geoms["geom"]
    circle = CircleSkyRegion(center=target_position, radius=1 * u.deg)
    regions = [circle]
    exclusion_mask = ~geom.region_mask(regions=regions)
    fov_bkg_maker = FoVBackgroundMaker(method="fit", exclusion_mask=exclusion_mask)

    gls = []
    bkg_norms = {}
    for dataset in datasets:
        bkg_model = FoVBackgroundModel(
            name=f"{dataset.name}_bkg", dataset_name=dataset.name
        )
        dataset.models = [bkg_model]
        dataset = fov_bkg_maker.run(dataset)
        bkg_norms[dataset.name] = bkg_model.parameters["norm"].value
        gl = GammapyLike(dataset.name, frame="galactic")
        gl.set_datasets(dataset)
        gl.set_background_models(bkg_model)
        gl.set_sources("rxj1713")
        gls.append(gl)

    pl = Powerlaw()
    spat = Gaussian_on_sphere(
        lon0=target_position.transform_to("galactic").l.deg,
        lat0=target_position.transform_to("galactic").b.deg,
        sigma=0.25,
    )
    es = ExtendedSource(source_name="rxj1713", spectral_shape=pl, spatial_shape=spat)
    pl.index.value = -2
    pl.index.prior = Uniform_prior(lower_bound=-3, upper_bound=-1)
    pl.K = 2.3 * 1e-11
    pl.K.prior = Log_uniform_prior(lower_bound=1e-18, upper_bound=1e-8)
    pl.piv.value = 1
    pl.piv.free = False
    spat.lon0.free = False
    spat.lat0.free = False
    spat.sigma.free = True
    spat.sigma.prior = Log_uniform_prior(lower_bound=0.1, upper_bound=1.0)
    model = Model(es)
    conv = AstromodelConverter(model, frame="galactic")
    for gl in gls:
        gl.set_model(model, conv)
        assert np.isclose(
            bkg_norms[gl.name],
            gl.nuisance_parameters[gl.name + "." + gl.name + "_bkg.norm"].value,
        )
