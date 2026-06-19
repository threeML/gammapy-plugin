---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: Python 3 
    language: python
---

# 3D Analysis of an Extended Source

In this example we fit the spectum and spatial extension of RX J1713.7-3946 (called rxj1713 in the following because typing is hard) first jointly and then stacked

## Joint Fit


Load all the relevant modules and set the energy unit to `u.TeV` as we are only dealing with H.E.S.S. data here

Additionally we load the `Gaussian_on_sphere` from `gammapy_plugin.utils.astromodels_functions` and `Cutoff_powerlaw` from `astromodels` for modelling the source later

```python
%matplotlib inline
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from astromodels.core.model import Model
from astromodels.core.units import get_units
from astromodels.functions import (
    Log_uniform_prior,
    Cutoff_powerlaw,
    Uniform_prior,
)
from astromodels.sources.extended_source import ExtendedSource
from astropy.coordinates import SkyCoord
from gammapy.data import DataStore
from gammapy.datasets import MapDataset, Datasets
from gammapy.makers import (
    MapDatasetMaker,
    SafeMaskMaker,
)
from gammapy.maps import MapAxis, WcsGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    PowerLawSpectralModel,
    SkyModel,
    GaussianSpatialModel,
    FoVBackgroundModel,
)
from regions import CircleSkyRegion
from threeML import DataList, load_analysis_results
from threeML.bayesian.bayesian_analysis import BayesianAnalysis
from threeML.classicMLE.joint_likelihood import JointLikelihood
from threeML.utils.progress_bar import trange

from gammapy_plugin.converter import AstromodelConverter
from gammapy_plugin.gammapy_like import GammapyLike
from astromodels.functions.functions_2D import Gaussian_on_sphere

get_units().energy = u.TeV
```

Select all the observations within a `5 deg` radius around RX J1713.7-3946

Prepare the geometry by setting the energy axis and `WcsGeom`

```python
datastore = DataStore.from_dir("$GAMMAPY_DATA/hess-dl3-dr1/")
```

```python
datastore.obs_table[datastore.obs_table["OBJECT"] == "MSH15-52"]
```

```python
target_position = SkyCoord.from_name("MSH15-52").galactic

selection = dict(
    type="sky_circle",
    frame="galactic",
    lon=target_position.l,
    lat=target_position.b,
    radius="5deg",
)
select_obs_tab = datastore.obs_table.select_observations(selection)

obs = datastore.get_observations(select_obs_tab["OBS_ID"])
print(len(obs))
# Prepare the geometry
energy_axis = MapAxis.from_energy_bounds(0.2, 40.0, 10, per_decade=True, unit="TeV")
energy_axis_true = MapAxis.from_energy_bounds(
    0.1, 100, 20, per_decade=True, unit="TeV", name="energy_true"
)
geom = WcsGeom.create(
    skydir=target_position,
    binsz=0.02,
    width=(6 * u.deg, 6 * u.deg),
    frame="galactic",
    axes=[energy_axis],
)
```

```python
geom
```

Create the relevant `Makers`

We will also run the `FoVBackgroundMaker` when jointly fitting - this is not needed

We exclude a `1 deg` circle at the source poition for fitting the background models to reduce overestimation of the background


```python
circle = CircleSkyRegion(center=target_position, radius=1 * u.deg)
circle1 = CircleSkyRegion(center=SkyCoord.from_name("HESS J1503-582"), radius=1 * u.deg)
circle2 = CircleSkyRegion(center=SkyCoord.from_name("HESS J1458-608"), radius=1 * u.deg)
circle3 = CircleSkyRegion(center=SkyCoord.from_name("HESS J1458-608"), radius=1 * u.deg)
regions = [circle, circle1, circle2, circle3]
exclusion_mask = ~geom.region_mask(regions=regions)
maker = MapDatasetMaker(
    selection=["counts", "background", "psf", "edisp", "exposure"],
)
safe_mask_maker = SafeMaskMaker(
    methods=["offset-max", "aeff-max", "bkg-peak"], offset_max="2.3 deg"
)
fov_bkg_maker = FoVBackgroundMaker(method="fit", exclusion_mask=exclusion_mask)
```

And now the datasets as well as the `GammapyPlugin` instances:

Few noteworthy things:
- we set the `FoVBackgroundModel` before and provide a custom `name`. As mentioned before, the background model is fitted afterwards and running the `FoVBackgroundMaker` is not necessary, but here we also print the fitted values from the `FoVBackgroundMaker` to compare to our results later
- in case you want to fit the background before and fix it during the fit you will always have to set the model before providing a non-default name as it defaults to `dataset_name + "-bkg"` which is not compatible with `astromodels`
- the `exclusion_mask` is only used for the `FoVBackgroundMaker` and will _not_ be used later during the acutal sampling
- the used map size of `6deg x 6deg` is fairly large and there are other sources in this region which we do not mask for simplicity: **_This will lead to overestimation of the background and therefore an underestimation of the source's flux!_**
- we set the dataset using default `mode="individual"` as we use on plugin per instance, afterwards we add the corresponding background model, which then gets parsed by `gammapy_plugin` and the normalization is treated as a nuisance parameter
- we used a galactic `WcsGeom` frame and we therefore also use galactic in the plugin

```python
datasets = Datasets()
gls = []
for o in obs:
    dataset = MapDataset.create(
        geom=geom, energy_axis_true=energy_axis_true, name=f"HESS_{o.obs_id}"
    )
    dataset = maker.run(dataset, o)
    dataset = safe_mask_maker.run(dataset, o)
    bkg_model = FoVBackgroundModel(name=f"{o.obs_id}_bkg", dataset_name=dataset.name)
    dataset.models = [bkg_model]
    dataset = fov_bkg_maker.run(dataset)
    print(
        f"Bkg norm for HESS_{o.obs_id}: {round(bkg_model.parameters['norm'].value,3)} +/- {round(bkg_model.parameters['norm'].error,3)}"
    )
    datasets.append(dataset)

gl = GammapyLike(dataset.name, frame="galactic")
gl.set_datasets(datasets, mode="stacked")
# gl.set_background_models(bkg_model)
gls.append(gl)
```

```python

```

## Setting up the Model

Using a exponential cutoff-powerlaw spectral model and a 2D Gaussian spatial one for simplicity

```python
cpl = Cutoff_powerlaw()
spat = Gaussian_on_sphere(
    lon0=target_position.transform_to("galactic").l.deg,
    lat0=target_position.transform_to("galactic").b.deg,
    sigma=0.5,
)
es = ExtendedSource(source_name="rxj1713", spectral_shape=cpl, spatial_shape=spat)
cpl.K = 5e-15 * u.Unit("TeV-1 cm-2 s-1")
cpl.index = -2 * u.dimensionless_unscaled
cpl.xc = 15 * u.TeV
cpl.piv = 1 * u.TeV
cpl.piv.free = False
cpl.K.prior = Log_uniform_prior(lower_bound=1e-15, upper_bound=1e-14)
cpl.index.prior = Uniform_prior(lower_bound=-3, upper_bound=-1)
cpl.xc.prior = Uniform_prior(lower_bound=1, upper_bound=30)
spat.lon0.free = False
spat.lat0.free = False
spat.sigma.unit = u.deg
spat.sigma.free = True
spat.sigma.prior = Uniform_prior(lower_bound=0.1, upper_bound=1)

model = Model(es)
```

We can either convert the model before setting it or a `AstromodelConverter` will be created when adding the model to plugins

```python
conv = AstromodelConverter(model, frame="galactic")
```

Now adding the model and the converted one to the plugins and specificly including the source. All sources are used by default, here we set it explicitly for no specific reason

```python
for gl in gls:
    gl.set_sources("rxj1713")
    gl.set_model(model, conv)
```

This next cell might take hours or days to finish using only a single core. Just skip it and load the result provided in the `result.fits` in the next cell.

```python
ba = BayesianAnalysis(model, DataList(*gls))
ba.set_sampler("multinest")
ba.sampler.setup()

res = ba.sample(quiet=False)
```

```python
# res = load_analysis_results("data/result.fits")
res.display()
```

Take a look at the corner plot of the source's free paramters:

```python
fig = res.corner_plot(
    components=[
        "rxj1713.Gaussian_on_sphere.sigma",
        "rxj1713.spectrum.main.Cutoff_powerlaw.K",
        "rxj1713.spectrum.main.Cutoff_powerlaw.index",
        "rxj1713.spectrum.main.Cutoff_powerlaw.xc",
    ]
)
fig.show()
```

Lets take a look at the differential flux at 1 TeV by integrating over the spatial part - the units of the normalzation `K` are actually `cm-2 TeV-1 s-1 deg-2` and not as stated in the results `cm-2 TeV-1 s-1`

```python
flux_1tev = res.optimized_model.extended_sources["rxj1713"].spectrum.main.shape(
    1 * u.TeV
) * res.optimized_model.extended_sources[
    "rxj1713"
].spatial_shape.get_total_spatial_integral(
    z=1 * u.TeV
)


samples = res.samples
n_samples = 1000
vals = np.zeros(n_samples)
cpl_copy = Cutoff_powerlaw()
spat_copy = Gaussian_on_sphere(
    lon0=target_position.transform_to("galactic").l.deg,
    lat0=target_position.transform_to("galactic").b.deg,
    sigma=0.5,
)
es_copy = ExtendedSource(
    source_name="dummy", spectral_shape=cpl_copy, spatial_shape=spat_copy
)

for i in trange(n_samples):
    n = np.random.choice(np.arange(samples.shape[-1]))
    spat_copy.sigma.value = samples[0, n]
    cpl_copy.K.value = samples[1, n]
    cpl_copy.index.value = samples[2, n]
    cpl_copy.xc.value = samples[3, n]
    vals[i] = cpl_copy(1) * spat_copy.get_total_spatial_integral()
hpd_1tev = np.array((np.percentile(vals, 2.5), np.percentile(vals, 97.5)))
print(f"Flux: {flux_1tev} +/- {np.abs((hpd_1tev-flux_1tev.value)[::-1])}")
print(f"95% hpd interval: {hpd_1tev}")
```

Keeping in mind that we are overestimating the background for sure, this is on the order of the value from [H.E.S.S. Collaboration, 2018](http://dx.doi.org/10.1051/0004-6361/201629790), who determined a flux of $2.3\pm0.1\:\mathrm{e}{-11}$ cm-2 s-1 TeV-1



```python
plt.matshow(res.get_correlation_matrix())
ticks = [
    i.replace("rxj1713.", "")
    .replace("_bkg_norm", "")
    .replace("spectrum.main.Cutoff_powerlaw.", "")
    .replace("HESS_", "")
    .replace("Gaussian_on_sphere.", "")
    .split("_")[0]
    for i in res.optimized_model.free_parameters.keys()
]
plt.xticks(np.arange(19), ticks, rotation=90)
plt.yticks(np.arange(19), ticks)
plt.colorbar()
```

```python
datasets_temp = Datasets()
for gl in gls:
    datasets_temp.append(gl.datasets[0])
```

Taking a look at the energy integrated counts map of these 15 observations, a 2D-Gaussian might not be the best description and a Powerlaw neither ;)

```python
datasets.stack_reduce().counts.sum_over_axes(keepdims=False).smooth(0.02 * u.deg).plot()
```

## Stacking the Dataset

We can also speed things a lot by
1. stacking the datasets
2. only fitting the FoVBackgroundModels before hand and fixing their values


```python
datasets_stacked = Datasets()

for o in obs:
    dataset = MapDataset.create(
        geom=geom, energy_axis_true=energy_axis_true, name=f"HESS_{o.obs_id}"
    )
    dataset = maker.run(dataset, o)
    dataset = safe_mask_maker.run(dataset, o)
    bkg_model = FoVBackgroundModel(name=f"{o.obs_id}_bkg", dataset_name=dataset.name)
    dataset.models = [bkg_model]
    dataset = fov_bkg_maker.run(dataset)
    print(
        f"Bkg norm for HESS_{o.obs_id}: {round(bkg_model.parameters['norm'].value,3)} +/- {round(bkg_model.parameters['norm'].error,3)}"
    )
    datasets_stacked.append(dataset)
```

```python jupyter={"outputs_hidden": true}
for d in datasets_stacked:
    print(d)
```

```python
datasets_stacked.stack_reduce().counts.sum_over_axes(keepdims=False).smooth(
    0.05 * u.deg
).plot()
```

Just stack them before passing them to the `GammapyLike` instance or set the `mode="stacked"` so the plugin does it for you :)

```python
gl_stacked = GammapyLike("stacked", frame="galactic")
gl_stacked.set_datasets(datasets_stacked, mode="stacked")
```

```python

```

Just the same model as before

```python
from astromodels.functions import Disk_on_sphere
```

```python
from astromodels import *
```

```python
cpl_stacked = Powerlaw()
spat_stacked = Gaussian_on_sphere(
    lon0=target_position.transform_to("galactic").l.deg,
    lat0=target_position.transform_to("galactic").b.deg,
    sigma=0.15,
)
es_stacked = ExtendedSource(
    source_name="rxj1713_stacked",
    spectral_shape=cpl_stacked,
    spatial_shape=spat_stacked,
)
cpl_stacked.K = 2.58e-12 / 3283 * u.Unit("TeV-1 cm-2 s-1")
cpl_stacked.index = -2.26 * u.dimensionless_unscaled
# cpl_stacked.xc = 8.5 * u.TeV
cpl_stacked.piv = 1 * u.TeV
cpl_stacked.piv.free = False
cpl_stacked.K.prior = Log_uniform_prior(lower_bound=1e-15, upper_bound=1e-8)
cpl_stacked.index.prior = Uniform_prior(lower_bound=-3, upper_bound=-1)
# cpl_stacked.xc.prior = Uniform_prior(lower_bound=1, upper_bound=30)
spat_stacked.lon0.free = False
spat_stacked.lat0.free = False
spat_stacked.sigma.free = True
spat_stacked.sigma.unit = u.deg
spat_stacked.sigma.prior = Log_uniform_prior(lower_bound=0.01, upper_bound=1)

model_stacked = Model(es_stacked)
conv_stacked = AstromodelConverter(model_stacked, "galactic")
```

```python
es_stacked
```

```python
gl_stacked.set_model(model_stacked, conv_stacked)
```

For simplicity just use minuit to minimize the likelihood

```python
gl_stacked.get_log_like()
```

```python
1509216 / 755
```

```python
from scipy.integrate import nquad

nquad(
    spat_stacked.evaluate,
    [(310, 330), (-10, 10)],
    args=(320.278110812368, -1.0728068, 0.15),
)
```

```python
(180 / np.pi) ** 2
```

```python
es_stacked.spectrum.main(1 * u.TeV) * es_stacked.spatial_shape(
    320.278, -1.072, 320.278, -1.072, 0.145
)
```

```python
gl_stacked.datasets[0]
```

```python
ba2 = BayesianAnalysis(model_stacked, DataList(gl_stacked))
```

```python
import time
```

```python

```

```python
ba2.set_sampler("multinest")
ba2.sampler.setup(verbose=True, resume=False)
fit = ba2.sample()
```

```python
res = ba2.results
```

```python
ba2.results.plot_chains()
```

```python
res.display()
```

```python
res.optimized_model.extended_sources["rxj1713_stacked"].get_spatially_integrated_flux(
    1 * u.TeV
)
```

```python
es_stacked.get_spatially_integrated_flux(1), es_stacked.get_spatially_integrated_flux(
    1 * u.TeV
)
```

```python
import time

start = time.time()
jl.fit()
stop = time.time()
print(f"Fit finished in {round(stop-start,3)}s")
res_stacked = jl.results
```

```python
res_stacked.display()
```

```python

```

```python
flux_1tev_stacked = res_stacked.optimized_model.extended_sources[
    "rxj1713_stacked"
].spectrum.main.shape(1 * u.TeV) * res_stacked.optimized_model.extended_sources[
    "rxj1713_stacked"
].spatial_shape.get_total_spatial_integral(
    z=1 * u.TeV
)

print(f"Flux: {flux_1tev_stacked}")
```

```python
es_stacked.spectrum.main(1 * u.TeV) * es_stacked.spatial_shape(
    target_position.l.deg, target_position.b.deg
)
```

```python
es_stacked.components["main"](1 * u.TeV) * (
    180 / 4 * np.pi
) ** 2  # /(0.5**2*np.pi)#*42000#*4*np.pi/(np.pi*0.5/180*np.pi)
```

```python
E = 1 * u.TeV
es_stacked.spatial_shape.get_total_spatial_integral(E) * es_stacked.spectrum.main.shape(
    E
)
```

```python

```

```python
spat_stacked.ndim
```

```python
es_stacked.get_spatially_integrated_flux(1 * u.TeV)
```

```python

```

```python

```

```python
e = np.geomspace(0.1, 50, 20) * u.TeV
plt.plot(e, es_stacked.get_spatially_integrated_flux(e))
plt.xscale("log")
plt.yscale("log")
```

```python
np.power(0.078, -1)
```

```python
es_stacked
```

```python
from astromodels import *
```

```python
gauss = Gaussian_on_sphere()
cpl = Cutoff_powerlaw()
es = ExtendedSource("test", spatial_shape=gauss, spectral_shape=cpl)
```

```python
cpl
```

```python
es.spectrum.main(2 * u.TeV)
```

```python
es.get_spatially_integrated_flux(2 * u.TeV)
```

```python

```

```python
E = 1 * u.TeV
es.spatial_shape.get_total_spatial_integral(E) * es.spectrum.main.shape(E)
```

Compared to the multiple hours `multinest` took to sample the full posterior using 10 threads this is way way faster at the cost of all the disadvantages of stacking the data ;)

```python
os.environ["GAMMAPY_DATA"]
```

```python

```

```python
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

# %matplotlib inline
from IPython.display import display
from gammapy.data import FixedPointingInfo, Observation, observatory_locations
from gammapy.datasets import MapDataset
from gammapy.irf import load_irf_dict_from_file
from gammapy.makers import MapDatasetMaker, SafeMaskMaker
from gammapy.maps import MapAxis, WcsGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    FoVBackgroundModel,
    GaussianSpatialModel,
    Models,
    PowerLawSpectralModel,
    SkyModel,
)

# Loading IRFs
irfs = load_irf_dict_from_file(
    "$GAMMAPY_DATA/cta-1dc/caldb/data/cta/1dc/bcf/South_z20_50h/irf_file.fits"
)
# Define the observation parameters (typically the observation duration and the pointing position):
livetime = 20.0 * u.hr
pointing_position = SkyCoord(ra=5, dec=0, unit="deg", frame="icrs")
# We want to simulate an observation pointing at a fixed position in the sky.
# For this, we use the `FixedPointingInfo` class
pointing = FixedPointingInfo(
    fixed_icrs=pointing_position.icrs,
)

# Define map geometry for binned simulation
energy_reco = MapAxis.from_edges(
    np.geomspace(0.5, 10, 20), unit="TeV", name="energy", interp="log"
)
geom = WcsGeom.create(
    skydir=pointing_position,
    binsz=0.02,
    width=(5, 5),
    frame="icrs",
    axes=[energy_reco],
)
# It is usually useful to have a separate binning for the true energy axis
energy_true = MapAxis.from_edges(
    np.geomspace(0.2, 20, 40), unit="TeV", name="energy_true", interp="log"
)

empty = MapDataset.create(geom, name="dataset-simu", energy_axis_true=energy_true)

# Define sky model to used simulate the data.
# Here we use a Gaussian spatial model and a Power Law spectral model.
spatial_model = GaussianSpatialModel(
    lon_0="4 deg", lat_0="0.1 deg", sigma="0.3 deg", frame="icrs"
)
spectral_model = PowerLawSpectralModel(
    index=3, amplitude="1e-11 cm-2 s-1 TeV-1", reference="1 TeV"
)
model_simu = SkyModel(
    spatial_model=spatial_model,
    spectral_model=spectral_model,
    name="model-simu",
)

bkg_model = FoVBackgroundModel(dataset_name="dataset-simu")

models = Models([model_simu, bkg_model])
# print(models)


######################################################################
# Now, comes the main part of dataset simulation. We create an in-memory
# observation and an empty dataset. We then predict the number of counts
# for the given model, and Poisson fluctuate it using ``fake()`` to make
# a simulated counts maps. Keep in mind that it is important to specify
# the ``selection`` of the maps that you want to produce
#

# Create an in-memory observation
location = observatory_locations["ctao_south"]
obs = Observation.create(
    pointing=pointing, livetime=livetime, irfs=irfs, location=location
)
# print(obs)

# Make the MapDataset
maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])

maker_safe_mask = SafeMaskMaker(methods=["offset-max"], offset_max=4.0 * u.deg)

dataset = maker.run(empty, obs)
dataset = maker_safe_mask.run(dataset, obs)
# print(dataset)

# Add the model on the dataset and Poisson fluctuate
dataset.models = models
dataset.fake()

fake_dataset = dataset.copy()
# Do a print on the dataset - there is now a counts maps
print(dataset)


dataset.counts.smooth(0.05 * u.deg).plot_interactive(add_cbar=True, stretch="linear")
plt.show()

models_fit = models.copy()

# We do not want to fit the background in this case, so we will freeze the parameters
models_fit["dataset-simu-bkg"].spectral_model.norm.frozen = True
models_fit["dataset-simu-bkg"].spectral_model.tilt.frozen = True

dataset.models = models_fit
print(dataset.models)

fit = Fit(optimize_opts={"print_level": 1})
result = fit.run(datasets=[dataset])

dataset.plot_residuals_spatial(method="diff/sqrt(model)", vmin=-0.5, vmax=0.5)
plt.show()


######################################################################
# Compare the injected and fitted models:
#

print(
    "True model: \n",
    model_simu,
    "\n\n Fitted model: \n",
    models_fit["model-simu"],
)


######################################################################
# Get the errors on the fitted parameters from the parameter table
#

display(result.parameters.to_table())
```

```python
gl_fake = GammapyLike(name="fake_ds", frame="icrs")
fake_dataset
fake_dataset.models = []
```

```python
gl_fake.set_datasets(fake_dataset)
```

```python
bkg_model_fake = FoVBackgroundModel(
    name="test_fake_bkg", dataset_name=fake_dataset.name
)
```

```python
gl_fake.set_background_models(bkg_model_fake)
```

```python
from astromodels import *
```

```python
pl = Powerlaw()
gauss = Gaussian_on_sphere()
es_fake = ExtendedSource(
    "fake_it_til_you_make_it", spatial_shape=gauss, spectral_shape=pl
)
```

```python
pl.K.min_value = 1e-20 * u.Unit("cm-2 s-1 TeV-1")
pl.K.max_value = 1e-5 * u.Unit("cm-2 s-1 TeV-1")
pl.K = 2e-11 * u.Unit("cm-2 s-1 TeV-1")
pl.K.set_uninformative_prior(Log_uniform_prior)
```

```python
pl.piv = 1 * u.TeV

pl.index.min_value = -5
pl.index.max_value = -1
pl.index.value = -2.5
pl.index.prior = Uniform_prior(lower_bound=-5, upper_bound=-1)
```

```python
gauss.lon0.min_value = 3.5
gauss.lon0.max_value = 4.5
gauss.lon0 = 4.1
gauss.lon0.free = True

gauss.lon0.set_uninformative_prior(Uniform_prior)


gauss.lat0.min_value = -1
gauss.lat0.max_value = 1
gauss.lat0 = 0
gauss.lat0.free = True
gauss.lat0.set_uninformative_prior(Uniform_prior)
gauss.sigma = 0.1 * u.deg
gauss.sigma.free = True
gauss.sigma.max_value = 2
gauss.sigma.prior = Log_uniform_prior(lower_bound=1e-2, upper_bound=1)
```

```python
model = Model(es_fake)
```

```python
model
```

```python
conv = AstromodelConverter(model, frame="icrs")
```

```python
gl_fake.set_model(model, conv)
```

```python
model
```

```python
model.fake_it_til_you_make_it.spectrum.main.shape.K.value = 9e-11
```

```python
print(gl_fake.get_log_like())
```

```python
lls = []
x = np.geomspace(0.01, 10000, 100)
for K in x:
    model.fake_it_til_you_make_it.spectrum.main.shape.K.value = K * 1e-11
    lls.append(gl_fake.get_log_like())
```

```python
plt.plot(x, lls)
plt.xscale("log")
```

```python
SkyCoord(ra=30, dec=0, unit="deg", frame="icrs").transform_to("galactic")
```

```python
# model.fake_it_til_you_make_it.spatial_shape.lat0.value = 0
```

```python
print(gl_fake.get_log_like())
gl_fake.gammapy_model["fake_it_til_you_make_it"]
```

```python
ba = BayesianAnalysis(model, data_list=DataList(gl_fake))
```

```python
ba.set_sampler("multinest")
ba.sampler.setup(verbose=True, resume=False)
```

```python
ba.sample()
```

```python
ba.results.corner_plot()
```

```python
gl_fake.datasets[0].counts.smooth(0.05 * u.deg).plot_interactive(
    add_cbar=True, stretch="linear"
)
```

```python
gl_fake.datasets[0].plot_residuals_spatial(
    method="diff/sqrt(model)", vmin=-0.5, vmax=0.5
)
```

```python
ba.results.plot_chains()
```

```python
import astropy.units as u
import numpy as np
from astromodels.core.model import Model
from astromodels.core.units import get_units
from astromodels.functions import (
    Gaussian_on_sphere,
    Log_uniform_prior,
    Powerlaw,
    Uniform_prior,
)
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
from gammapy_plugin.gammapy_plugin import GammapyLike
from gammapy_plugin.test.utils import read_in_gammapy_datasets
from gammapy_plugin.utils.package_data import get_path_of_data_dir

get_units().energy = u.TeV
target_position = SkyCoord.from_name("RX J1713.7-3946").galactic


datasets = read_in_gammapy_datasets(
    get_path_of_data_dir().joinpath("test/rxj17137_3946/")
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
```


```python
jl = JointLikelihood(model, DataList(gl))
jl.fit()
res = jl.results
```


```python

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
```


```python
es
```

```python
E = np.array([1]) * u.TeV
```

```python
es(
    np.array([347.26942958693814]) * u.deg,
    np.array([-0.2569129018602152]) * u.deg,
    np.array([1, 1]) * u.TeV,
).to(u.Unit("TeV-1 s-1 cm-2 deg-2"))
```

```python
resu.models
```

```python

```
