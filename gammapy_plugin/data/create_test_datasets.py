import logging
import os
import astropy.units as u
from astropy.coordinates import SkyCoord
from gammapy.data import DataStore
from gammapy.datasets import Datasets, MapDataset
from gammapy.makers import FoVBackgroundMaker, MapDatasetMaker, SafeMaskMaker
from gammapy.maps import MapAxis, WcsGeom
from regions import CircleSkyRegion

from gammapy_plugin.utils.package_data import get_path_of_data_dir

log = logging.getLogger(__name__)

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
for o in obs:
    dataset = MapDataset.create(
        geom=geom, energy_axis_true=energy_axis_true, name=f"HESS_{o.obs_id}"
    )
    dataset = maker.run(dataset, o)
    dataset = safe_mask_maker.run(dataset, o)
    datasets.append(dataset)
base_dir = get_path_of_data_dir().joinpath("test/rxj17137_3946/")
if not os.path.exists(base_dir):
    os.makedirs(base_dir)
for d in datasets:
    fn = str(base_dir.joinpath(f"{d.name}.fits"))
    log.warning(f"Saving to {fn}")
    d.write(
        filename=fn,
        overwrite=True,
    )
