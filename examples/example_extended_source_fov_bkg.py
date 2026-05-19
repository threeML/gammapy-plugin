import astropy.units as u
from astromodels.core.model import Model
from astromodels.core.units import get_units
from astromodels.functions import Log_uniform_prior, Powerlaw, Uniform_prior
from astromodels.sources.extended_source import ExtendedSource
from astropy.coordinates import SkyCoord
from gammapy.data import DataStore
from gammapy.datasets import Datasets, MapDataset
from gammapy.makers import FoVBackgroundMaker, MapDatasetMaker, SafeMaskMaker
from gammapy.maps import MapAxis, WcsGeom
from gammapy.modeling.models import FoVBackgroundModel
from mpi4py.MPI import COMM_WORLD
from regions import CircleSkyRegion
from threeML import BayesianAnalysis, DataList

from gammapy_plugin.converter import AstromodelConverter
from gammapy_plugin.gammapy_like import GammapyLike
from gammapy_plugin.utils.astromodels_functions import Gaussian_on_sphere

comm = COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

get_units().energy = u.TeV


datastore = DataStore.from_dir("$GAMMAPY_DATA/hess-dl3-dr1/")
target_position = SkyCoord.from_name("RX J1713.7-3946").galactic

selection = dict(
    type="sky_circle",
    frame="galactic",
    lon=target_position.l,
    lat=target_position.b,
    radius="5deg",
)
select_obs_tab = datastore.obs_table.select_observations(selection)

obs = datastore.get_observations(select_obs_tab["OBS_ID"])

# Prepare the geometry
energy_axis = MapAxis.from_energy_bounds(0.3, 10.0, 15, unit="TeV")
energy_axis_true = MapAxis.from_energy_bounds(
    0.1, 20, 10, per_decade=True, unit="TeV", name="energy_true"
)
geom = WcsGeom.create(
    skydir=target_position,
    binsz=0.02,
    width=(6 * u.deg, 6 * u.deg),
    frame="galactic",
    axes=[energy_axis],
)

circle = CircleSkyRegion(center=target_position, radius=1 * u.deg)
regions = [circle]
exclusion_mask = ~geom.region_mask(regions=regions)
maker = MapDatasetMaker(
    selection=["counts", "background", "psf", "edisp", "exposure"],
)
safe_mask_maker = SafeMaskMaker(
    methods=["offset-max", "aeff-max", "bkg-peak"], offset_max="2.3 deg"
)
fov_bkg_maker = FoVBackgroundMaker(method="fit", exclusion_mask=exclusion_mask)


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
    datasets.append(dataset)
    gl = GammapyLike(dataset.name, frame="galactic")
    gl.set_datasets(dataset)
    gl.set_background_models(bkg_model)
    gls.append(gl)

pl = Powerlaw()
spat = Gaussian_on_sphere(
    lon0=target_position.transform_to("galactic").l.deg,
    lat0=target_position.transform_to("galactic").b.deg,
    sigma=0.5,
)
es = ExtendedSource(source_name="rxj1713", spectral_shape=pl, spatial_shape=spat)
pl.index.value = -2
pl.index.prior = Uniform_prior(lower_bound=-2.5, upper_bound=-1.9)
pl.K = 8 * 1e-15
pl.K.prior = Log_uniform_prior(lower_bound=1e-16, upper_bound=1e-14)
pl.piv.value = 1
pl.piv.free = False
spat.lon0.free = False
spat.lat0.free = False
spat.sigma.free = True
spat.sigma.prior = Uniform_prior(lower_bound=0.25, upper_bound=0.75)
model = Model(es)

conv = AstromodelConverter(model, frame="galactic")
for gl in gls:
    gl.set_sources("rxj1713")
    gl.set_model(model, conv)

comm.Barrier()
ba = BayesianAnalysis(model, DataList(*gls))
ba.set_sampler("ultranest")
ba.sampler.setup()


ba.sample(quiet=False)
if rank == 0:
    res = ba.results
    res.write_to("result_go_big.fits", overwrite=True)
