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
Model to a gammapy SkyModel, so we can perform the forward folding of said model in 
gammapy.
Let's have a look!

We start by importing everything for a simple astromodels Model
```python
from astromodels.core.model import Model
from astromodels.functions import Powerlaw
from astromodels.functions import Disk_on_sphere
from astromodels.sources import ExtendedSource,PointSource
from astromodels.utils.io import display
import astropy.units as u
```

which we then define
```python
pl = Powerlaw()
disk = Disk_on_sphere()
source = ExtendedSource("test_source",spectral_shape=pl,spatial_shape=disk)
disk.lon0 = 1*u.deg
disk.lat0 = 10*u.deg
model = Model(source)
```

and take a look at what we created
```python
display(model)
```

B.E.A.U.T.I.F.U.L.

Ok, now for the conversion - first we need to import stuff again:

```python
from gammapy_plugin.converter import AstromodelConverter
```

and it is indeed as simple as running

```python
conv = AstromodelConverter(model = model)
```

We now have an instance of the converter which has already perfomed the conversion 
up on initialization.
We only created one single source so we only have one single model which we can index

```python
conv.gammapy_models[0]
```

Nice!
We see that both our `lon0`and `lat0` values are the ones we assigned before :)

Now let's update them. If we simply run

```python
pl.K.value = 10
conv.gammapy_models[0]
```

We see that nothing has changed :/
We have to invoke the `update()` function:
```python
conv.update()
conv.gammapy_models[0]
```

Perfect!

This calls the indiviudal update methods of the sources which then update all the models.

## Point Sources

In case you are only analysing a spectrum dataset and do not care about pixelatiion on a
map or the PSF at all you can set the `convert_ps` flag to `False`.

By default this is set to true, meaning when you transform an astromodels `PointSource`
you will get a `PointSpatialModel` (see e.g. 
[here](https://docs.gammapy.org/2.1/user-guide/model-gallery/spatial/plot_point.html)):

```python
pl = Powerlaw()
ps = PointSource("test_ps",spectral_shape=pl,ra = 0, dec =12.3 )
model_ps = Model(ps)
conv_ps = AstromodelConverter(model_ps)
conv_ps.gammapy_models[0]
```

If you set it to `False` you will not get any spatial component:

```python
pl = Powerlaw()
ps = PointSource("test_ps",spectral_shape=pl,ra = 0, dec = 12.3)
model_ps = Model(ps)
conv_ps = AstromodelConverter(model_ps,convert_ps = False)
conv_ps.gammapy_models[0]
```


## Plotting
Stealing from the 
[gammapy docs](https://docs.gammapy.org/2.1/user-guide/model-gallery/spectral/plot_powerlaw.html)
we can now also use the same plotting routines for the models:

But first let's update our parameters to something in the "`gammapy` energy range"
```python
pl.K = 1*u.Unit("TeV-1 cm-2 s-1")
pl.piv = 1*u.TeV
```
and please be aware that `astromodels` uses $K\times\left(\frac{E}{piv}\right)^{i}$
not $K\times\left(\frac{E}{piv}\right)^{-i}$ for a powerlaw definition.

```python
from astropy import units as u
import matplotlib.pyplot as plt


energy_bounds = [0.1, 100] * u.TeV
conv_ps.gammapy_models[0].spectral_model.plot(energy_bounds)
plt.grid(which="both")

```

We can also do the same for the Disk spatial model before
```python
from gammapy.maps import Map, WcsGeom

pl = Powerlaw()
disk = Disk_on_sphere()
source = ExtendedSource("test_source",spectral_shape=pl,spatial_shape=disk)
disk.lon0 = 20*u.deg
disk.lat0 = 0*u.deg
disk.radius = 1*u.deg
model = Model(source)
conv = AstromodelConverter(model = model)

lon_0 = disk.lon0.value
lat_0 = disk.lat0.value
reval = 2*disk.radius.value
dr = 0.02
geom = WcsGeom.create(
    skydir=(lon_0, lat_0),
    binsz=dr,
    width=(2 * reval, 2 * reval),
    frame="icrs",
)


fig, ax = plt.subplots(1, figsize=(9, 6))
meval = conv.gammapy_models[0].spatial_model.evaluate_geom(geom)
Map.from_geom(geom=geom, data=meval.value, unit=meval.unit).plot(ax=ax)
plt.tight_layout()
```
