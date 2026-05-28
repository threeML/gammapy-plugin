---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.17.3
  kernelspec:
    display_name: threeML
    language: python
    name: python3
    
---
# Spectrum Fitting 

In this first basic example we will perform a basic spectral analysis of the Crab nebula
using the public H.E.S.S. data release (already included in `gammapy`).

Let's get started by importing all the relevant stuff.

```python
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord

from astromodels.core.model import Model
from astromodels.core.units import get_units
from astromodels.functions import Log_parabola, Log_uniform_prior, Uniform_prior
from astromodels.sources import PointSource

from gammapy.data import DataStore
from gammapy.datasets import Datasets, SpectrumDataset
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import LogParabolaSpectralModel, SkyModel

from regions import CircleSkyRegion

from threeML import BayesianAnalysis
from threeML.data_list import DataList

from gammapy_plugin.converter import AstromodelConverter
from gammapy_plugin.gammapy_like import GammapyLike
from gammapy_plugin.test.utils import get_close
```

We will now set the `astromodels` energy unit to `TeV`.
This is a feature that requires `astromodels >= 2.5.1`

```python
get_units().energy = u.TeV
```

Alright perfect, let's continue by laoding the relevant data

```python
datastore = DataStore.from_dir("$GAMMAPY_DATA/hess-dl3-dr1/")
obs_ids = [23523, 23526, 23559, 23592]
observations = datastore.get_observations(obs_ids)
```

We perform a standard `gammapy`-workflow of creating a dataset

- setting as target
- excluding the target region
- creating our energy axis

```python
target_position = SkyCoord(ra=83.63, dec=22.01, unit="deg", frame="icrs")

on_region_radius = Angle("0.11 deg")
on_region = CircleSkyRegion(center=target_position.galactic, radius=on_region_radius)
exclusion_region = CircleSkyRegion(
    center=SkyCoord(183.604, -8.708, unit="deg", frame="galactic"),
    radius=0.5 * u.deg,
)

skydir = target_position.galactic
geom = WcsGeom.create(
    npix=(250, 250), binsz=0.02, skydir=skydir, proj="TAN", frame="galactic"
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
```

We now have everything to create our `makers` and run them.
In this example we use a `ReflectedRegionsBackgroundMaker`

```python
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
datasets_copy = datasets.copy()
```

Let's continue with the `threeML` steps
- choosing and setting up a model


```python
logp = Log_parabola()
ps = PointSource(
    source_name="crab",
    ra=target_position.ra.deg,
    dec=target_position.dec.deg,
    spectral_shape=logp,
)
logp.K.prior = Log_uniform_prior(lower_bound=1e-13, upper_bound=1e-9)
logp.K = 1e-11 * u.Unit("TeV-1 cm-2 s-1")
logp.piv = 1 * u.TeV
logp.piv.free = False
logp.alpha.prior = Uniform_prior(lower_bound=-3, upper_bound=-1)
logp.alpha = -2
logp.beta.prior = Uniform_prior(lower_bound=0, upper_bound=2)
logp.beta = 1

model = Model(ps)
```

```python
conv = AstromodelConverter(model, frame="galactic")
gl = GammapyLike("hess", sources="crab")
gl.set_datasets(datasets)
gl.set_model(model, converted_model=conv)

ba = BayesianAnalysis(model, DataList(gl))
ba.set_sampler("multinest")
ba.sampler.setup(resume = False)
ba.sample()
res = ba.results
res
```

```python tags=["nbsphinx-thumbnail"]
_ = res.corner_plot()
```
