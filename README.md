[![CI](https://github.com/threeML/gammapy-plugin/actions/workflows/build_test.yml/badge.svg)](https://github.com/threeML/gammapy-plugin/actions/workflows/build_test.yml)

# Gammapy Plugin
`threeML` plugin to use `gammapy` datasets in your `threeML` analysis.

This version is based on previous work by [J. Michael Burgess](https://github.com/grburgess), [Sajan Kumar](https://github.com/skumarudel) and [AnjanaTel](https://github.com/AnjanaTel).

## Installation
```bash
git clone https://github.com/threeML/gammapy-plugin.git
cd gammapy-plugin
pip install .
```

### Requirements

## Usage and Examples
For an example check out e.g. [the example notebook for an extened source](./examples/example_extended_source_fov_bkg.ipynb).

The basic procedure after creating a `gammapy` dataset is

```python
gl = GammapyLike(name = "name_of_the_plugin")   # initializing the plugin
gl.set_datasets(datasets,mode="individual")     # adding the gammapy dataset
gl.set_sources("name_of_the_source")            # setting the source
gl.set_model(model)                             # setting the astromodels model
```

When using `MapDataset` with  `FoVBackgroundModel` you can also include them in a fit.
The parameters of the background models are currently treated as nuisance parameters during the fitting.

This might then look like this for example ():

```python
datasets = Datasets()
gls = []
for o in obs:
    dataset = MapDataset.create(
        geom=geom, energy_axis_true=energy_axis_true, name=f"HESS_{o.obs_id}"
    )
    dataset = maker.run(dataset, o)
    dataset = safe_mask_maker.run(dataset, o)
    bkg_model = FoVBackgroundModel(name = f"{o.obs_id}_bkg",dataset_name= dataset.name)
    dataset.models = [bkg_model]
    dataset = fov_bkg_maker.run(dataset)
    datasets.append(dataset)
    gl = GammapyLike(dataset.name, frame="galactic")
    gl.set_datasets(dataset)
    gl.set_background_models(bkg_model)
    gls.append(gl)
```

For this [PR #5747](https://github.com/gammapy/gammapy/pull/5747) is needed additionally `gammapy==1.3`