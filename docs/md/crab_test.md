---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.7.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```python
from pathlib import Path
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion

from gammapy.data import DataStore
from gammapy.datasets import (
    Datasets,
    FluxPointsDataset,
    SpectrumDataset,
    SpectrumDatasetOnOff,
)
from gammapy.estimators import FluxPointsEstimator
from gammapy.estimators.utils import resample_energy_edges
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,
)
from gammapy.maps import MapAxis, RegionGeom, WcsGeom
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,
    SkyModel,
    create_crab_spectral_model,
)
```


```python
from gammapy_plugin.GammapyLike import GammapyLike
from threeML.bayesian.bayesian_analysis import BayesianAnalysis
from threeML.data_list import DataList
from astromodels.core.model import Model
from astromodels.sources import PointSource
from astromodels.functions import Log_parabola,Log_uniform_prior,Uniform_prior
from astromodels.functions.function import Function1D, FunctionMeta
```


```python
datastore = DataStore.from_dir("$GAMMAPY_DATA/hess-dl3-dr1/")
obs_ids = [23523, 23526, 23559, 23592]
observations = datastore.get_observations(obs_ids)
```


```python
target_position = SkyCoord(ra=83.63, dec=22.01, unit="deg", frame="icrs")
on_region_radius = Angle("0.11 deg")
on_region = CircleSkyRegion(center=target_position, radius=on_region_radius)
```


```python
exclusion_region = CircleSkyRegion(
    center=SkyCoord(183.604, -8.708, unit="deg", frame="galactic"),
    radius=0.5 * u.deg,
)

skydir = target_position.galactic
geom = WcsGeom.create(
    npix=(150, 150), binsz=0.05, skydir=skydir, proj="TAN", frame="icrs"
)

exclusion_mask = ~geom.region_mask([exclusion_region])

```


```python
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
```


```python
datasets = Datasets()

for obs_id, observation in zip(obs_ids, observations):
    dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
    dataset_on_off = bkg_maker.run(dataset, observation)
    dataset_on_off = safe_mask_maker.run(dataset_on_off, observation)
    datasets.append(dataset_on_off)


```


```python
Log_parabola?
```


```python
logp = Log_parabola()
logp.K.prior = Log_uniform_prior(lower_bound = 1e-23,upper_bound = 1e-15)
logp.piv.value = 1e9
logp.piv.free = False
logp.alpha.prior = Uniform_prior(lower_bound = -5,upper_bound = 5)
logp.beta.prior = Uniform_prior(lower_bound = 0,upper_bound = 2)
ps = PointSource(source_name = "crab",ra = target_position.ra.deg,dec = target_position.dec.deg,spectral_shape=logp)
```


```python
model = Model(ps)
```


```python
gl = GammapyLike("hess",sources = "crab")
gl.set_datasets(datasets)
gl.set_model(model)
```

```python
ba = BayesianAnalysis(model, data_list = DataList(gl))
ba.set_sampler("multinest")
ba.sampler.setup()
ba.sample()
```

```python
ba.results.plot_chains()
```


```python
ba.results.corner_plot()
```


```python
x = ba.results.optimized_model.crab.spectrum.main.Log_parabola.K
```


```python
ba.results.get_highest_density_posterior_interval("crab.spectrum.main.Log_parabola.K",cl=0.95)
```


```python
(x.value*x.unit).to("1/(TeV cm2 s)")
```


```python
from threeML.io.plotting.model_plot import plot_spectra
```


```python
y = plot_spectra(ba.results,ene_min = 0.1,ene_max = 100,energy_unit = "TeV",flux_unit = "TeV2/(TeV cm2 s)")
```


```python
y.gca().plot(np.geomspace(1e-1,1e2,100),np.power(np.geomspace(1e-1,1e2,100),2)*Log_parabola().evaluate(np.geomspace(1e-1,1e2,100),3.23*1e-11,1,-2.47,0.24),label = "MAGIC")
```


```python
y.gca().legend()
y
```


```python
from gammapy.modeling.models import LogParabolaSpectralModel
from gammapy.modeling import Fit
```


```python
logp_gammapy = LogParabolaSpectralModel(amplitude=1e-12*u.Unit("cm-2 s-1 TeV-1"),reference =1*u.TeV )
```


```python
models = SkyModel(spectral_model=logp_gammapy,name = "crab_gp")
dataset_stacked = Datasets(datasets).stack_reduce()
dataset_stacked.models = models
fit_stacked = Fit()
result_stacked = fit_stacked.run([dataset_stacked])
```


```python
models.spectral_model.to_dict()
```


```python
ba.results.estimate_covariance_matrix()
```


```python
ba.results.display()
```


```python
x = ba.results.get_median_fit_model().free_parameters["crab.spectrum.main.Log_parabola.K"]
```


```python
x
```


```python
np.isclose((x.value*x.unit).to("TeV-1 s-1 cm-2").value,models.spectral_model.to_dict()["spectral"]["parameters"][0]["value"])
```


```python
np.isclose(ba.results.get_median_fit_model().free_parameters["crab.spectrum.main.Log_parabola.beta"].value,models.spectral_model.to_dict()["spectral"]["parameters"][3]["value"])
```


```python
ba.results.get_data_frame(error_type = "hpd")
```


```python
models.spectral_model.to_dict()
```


```python
def get_close(threeML_results,gammapy_result_dict):
    try:
        bm = threeML_results.get_median_fit_model().free_parameters
    except AttributeError:
        bm = threeML_results.optimized_model.free_parameters
    hdp = threeML_results.get_data_frame(error_type = "hpd")
    for p in bm.keys():
        pn= p.split(".")[-1]
        if pn == "K":
            pn = "amplitude"
        for gp in gammapy_result_dict["spectral"]["parameters"]:
            if gp["name"] == pn:
                break
        if pn == "amplitude":
            v = (bm[p].value*bm[p].unit).to("TeV-1 cm-2 s-1").value
            min_v = ((hdp.loc[p]["negative_error"]+bm[p].value)*bm[p].unit).to("TeV-1 cm-2 s-1").value
            max_v = ((hdp.loc[p]["positive_error"]+bm[p].value)*bm[p].unit).to("TeV-1 cm-2 s-1").value
        elif pn == "alpha":
            v = -bm[p].value
            min_v = hdp.loc[p]["negative_error"]- bm[p].value
            max_v =hdp.loc[p]["positive_error"] - bm[p].value 
        else:
            v = bm[p].value
            min_v = hdp.loc[p]["negative_error"] +bm[p].value
            max_v = hdp.loc[p]["positive_error"] +bm[p].value
        return bool(gp["value"]<= max_v and min_v<=gp["value"])
```


```python
get_close(ba.results,models.spectral_model.to_dict())
```


```python
from threeML import JointLikelihood
```


```python
JointLikelihood?
```


```python
jl = JointLikelihood(model,DataList(gl))
```


```python
res = jl.fit()
```


```python
get_close(jl.results,models.spectral_model.to_dict())
```


```python
jl.results
```


```python
ba.results
```


```python
from importlib.metadata import version
version("gammapy")
```


```python

```
