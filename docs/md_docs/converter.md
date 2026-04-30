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
    name: threeml
    
---
# Converting a astromodels Model 

The core concept behind that plugin is the model converter that maps a astromodels
Model to a gammapy SkyModel.
Let's have a look:

We start by importing everything for a simple astromodels Model
```python
from astromodels.core.model import Model
from astromodels.functions import Powerlaw
from astromodels.functions import Disk_on_sphere
from astromodels.sources import ExtendedSource
```

which we then define
```python
pl = Powerlaw()
disk = Disk_on_sphere()
source = ExtendedSource("test_source",spectral_shape=pl,spatial_shape=disk)
model = Model(source)
```

