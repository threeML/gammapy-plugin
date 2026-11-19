---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.6
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Using astromodels in gammapy

In case you do not want to use `threeML` for fitting or sampling but still want to use
`astromodels` for modelling you can use the [converter](./converter.ipynb) to do so.

Let's start by defining an `astromodels` model for the Crab Nebula and doing some experimental `astromodels` stuff:


The default energy unit in `astromodels` is `keV` - which makes all values quiet small when having data in the `TeV` range.
There is currently no fully supported workaround for this and the easiest thing is to just be aware of this while modelling:
multiply your energies by `1e9` and your differential photon fluxes `1e-9` and you should be good.

Also consider setting the `min_value` and `max_value` to reasonably small values so the fitting engine will not have a hard time:

```python
%matplotlib inline
from astromodels.core.model import Model
from astromodels.functions import Powerlaw
from astromodels.sources import PointSource

pl = Powerlaw()
ps = PointSource("crab", ra=83.63, dec=22.01, spectral_shape=pl)
pl.piv.value = 1e9  # 1 TeV
pl.K.value = 1e-21  # 1e-12 TeV-1 cm-2 s-1
pl.K.min_value = 1e-23
pl.K.max_value = 1e-17
pl.index = -2  # remember the different powerlaw definitions
model = Model(ps)
```


Okay we are all set and can convert our model:

```python
from gammapy_plugin.converter import AstromodelConverter

conv = AstromodelConverter(model)
```

Let's do all the `gammapy` stuff - that is highly inspired by their [1D Spectral Analysis Tutorial](https://docs.gammapy.org/2.1/tutorials/analysis-1d/spectral_analysis.html) so you can compare the result
```python
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.coordinates import Angle, SkyCoord
from gammapy.data import DataStore
from gammapy.datasets import (
    Datasets,
    FluxPointsDataset,
    SpectrumDataset,
)
from gammapy.estimators import FluxPointsEstimator
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from gammapy.modeling import Fit
from IPython.display import display
from regions import CircleSkyRegion

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
    0.1, 40, nbin=10, per_decade=True, unit="TeV", name="energy"
)
energy_axis_true = MapAxis.from_energy_bounds(
    0.05, 100, nbin=20, per_decade=True, unit="TeV", name="energy_true"
)

geom = RegionGeom.create(region=on_region, axes=[energy_axis])
dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

dataset_maker = SpectrumDatasetMaker(
    containment_correction=True, selection=["counts", "exposure", "edisp"]
)
bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)
safe_mask_maker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=10)

datasets = Datasets()

for obs_id, observation in zip(obs_ids, observations):
    dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
    dataset_on_off = bkg_maker.run(dataset, observation)
    dataset_on_off = safe_mask_maker.run(dataset_on_off, observation)
    datasets.append(dataset_on_off)
print(datasets)
```
Now that we have our `Datasets` we can assign the `SkyModel`s we have created using the `AstromodelsConverter` and run the fit using
`gammapy` after that.

```python
datasets.models = conv.gammapy_models

fit_joint = Fit()
result_joint = fit_joint.run(datasets=datasets)
```

Moment of truth... did we converge?

```python
print(result_joint)
```

We did! &#x1F389;

And the parameter results also look good:

```python
display(result_joint.models.to_parameters_table())
```

We can also inspect our `gammapy` result as usual

```python
ax_spectrum, ax_residuals = datasets[0].plot_fit()
ax_spectrum.set_ylim(0.1, 40)
plt.show()
```
We can also calculate the flux points for this model

(please ignore all those "The current value of the parameter K (1.0) was above the new maximum 1e-17." warnings that may pop up - that's a `astromodels` issue - everything else works like it should)

```python
fpe = FluxPointsEstimator(
    energy_edges=energy_axis.edges, source="crab", selection_optional="all"
)
flux_points = fpe.run(datasets=datasets)
```

```python
fig, ax = plt.subplots()
flux_points.plot(ax=ax, sed_type="e2dnde", color="darkorange")
flux_points.plot_ts_profiles(ax=ax, sed_type="e2dnde")
ax.set_xlim(0.6, 40)
plt.show()
```

```python
from gammapy.modeling.models import Models
```

```python
flux_points_dataset = FluxPointsDataset(
    data=flux_points, models=Models(conv.gammapy_models[0].copy())
)
ax, _ = flux_points_dataset.plot_fit()
ax.set_xlim(0.6, 40)
plt.show()
```

So we can also model in `astromodels` and fit in `gammapy`.

The `gammapy-plugin` also is running tests in the CI on this, including `XSPEC` models via `astromodels`.
