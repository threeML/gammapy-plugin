[![CI](https://github.com/threeML/gammapy-plugin/actions/workflows/build_test_pip.yml/badge.svg?branch=main)](https://github.com/threeML/gammapy-plugin/actions/workflows/build_test_pip.yml)
[![Docs](https://app.readthedocs.org/projects/gammapy-plugin/badge/?version=latest&style=flat)](https://gammapy-plugin.readthedocs.io/en/latest/)
[![codecov](https://codecov.io/github/threeML/gammapy-plugin/graph/badge.svg?token=EGG9OUMT9J)](https://codecov.io/github/threeML/gammapy-plugin)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/threeML/gammapy-plugin/main.svg)](https://results.pre-commit.ci/latest/github/threeML/gammapy-plugin/main)

# Gammapy Plugin
![A gammapy plugin for threeML](docs/_static/banner-dark.svg)

`threeML` plugin to use `gammapy` datasets in your `threeML` analysis.


## Docs

You can find the docs containg also example notebooks 
either on [readthedocs](https://gammapy-plugin.readthedocs.io/en/latest/) or on the 
[github pages](https://threeml.github.io/gammapy-plugin/).

## Installation

The following should be sufficient and install all the relevant dependencies.

```bash
git clone https://github.com/threeML/gammapy-plugin.git
cd gammapy-plugin
pip install .
```

You can test the installation using `pytest`. To get all the dependencies it is
recommeneded to run 
```bash
# install test dependenices
pip install "gammapy_plugin[tests]"

# equivalent to gammapy download datasets - no need to do that if you already have them
# downloaded and stored in  $GAMMAPY_DATA
export GAMMAPY_DATA=/path/to/gammapy_datasets/
git clone https://github.com/gammapy/gammapy-data.git $GAMMAPY_DATA

# run pytests
python -m pytest --pyargs gammapy_plugin -n auto --dist worksteal 
```

The first part is equivalent to running `gammapy download datasets` [(see the `gammapy` installation
instructions)](https://docs.gammapy.org/1.3/getting-started/index.html#getting-started).

These data are required by the `pytest` tests.

### Requirements

Python 3.12 or 3.13 together with `gammapy>=v2.0.0`

## Usage and Examples

The basic procedure after creating a `gammapy` dataset is

```python
gl = GammapyLike(name="name_of_the_plugin")  # initializing the plugin
gl.set_datasets(datasets, mode="individual")  # adding the gammapy dataset
gl.set_sources("name_of_the_source")  # setting the source
gl.set_model(model)  # setting the astromodels model
```

When using `MapDataset` with `FoVBackgroundModel` you can also include them in a fit.
The parameters of the background models are currently treated as nuisance parameters during the fitting/sampling.

This might then look like this for example:

```python
datasets = Datasets()
gls = []
for o in observations:
    dataset = MapDataset.create(
        geom=geom, energy_axis_true=energy_axis_true, name=f"HESS_{o.obs_id}"
    )
    dataset = maker.run(dataset, o)
    dataset = safe_mask_maker.run(dataset, o)
    bkg_model = FoVBackgroundModel(name=f"{o.obs_id}_bkg", dataset_name=dataset.name)
    dataset.models = [bkg_model]
    dataset = fov_bkg_maker.run(dataset)
    datasets.append(dataset)
    gl = GammapyLike(dataset.name, frame="galactic")
    gl.set_datasets(dataset)
    gl.set_background_models(bkg_model)
    gls.append(gl)
```

# Notes
This current version is based on previous work by [J. Michael Burgess](https://github.com/grburgess), [Sajan
Kumar](https://github.com/skumarudel) and [AnjanaTel](https://github.com/AnjanaTel).

