import os
import glob
import pickle
import copy
import re
from astropy.io import fits
import matplotlib.pyplot as plt 
import sys, getopt
import yaml
import gammapy
import numpy as np
import astropy
import regions
from astropy.time import Time
#from pathlib import Path
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from regions import CircleSkyRegion
from gammapy.maps import Map, MapAxis, WcsGeom,  RegionGeom
from gammapy.modeling import Fit 
from gammapy.data import DataStore
from gammapy.datasets import (
    Datasets,
    SpectrumDataset,
    SpectrumDatasetOnOff,
    FluxPointsDataset,
    MapDataset,
)
from gammapy.modeling.models import (
    PowerLawSpectralModel,
    ExpCutoffPowerLawSpectralModel,
    #create_crab_spectral_model,
    SkyModel,
)
from gammapy.makers import (
    SafeMaskMaker,
    SpectrumDatasetMaker,
    ReflectedRegionsBackgroundMaker,
    MapDatasetMaker,
)


from threeML.exceptions.custom_exceptions import custom_warnings
from threeML.io.file_utils import sanitize_filename
from threeML.plugin_prototype import PluginPrototype
from threeML.utils.statistics.gammaln import logfactorial
from threeML.utils.unique_deterministic_tag import get_unique_deterministic_tag
from threeML.utils.power_of_two_utils import is_power_of_2
from threeML.io.package_data import get_path_of_data_file
from threeML.io.dict_with_pretty_print import DictWithPrettyPrint
from threeML.io.logging import setup_logger
from threeML.io.logging import setup_logger
log = setup_logger(__name__)

__instrument_name = "VERITAS (with gammapy)"



class GammapyLike(PluginPrototype):

    def __new__(cls, *args, **kwargs):
        instance = object.__new__(cls)
        return instance

    
    def __init__(self, name, config_file="config.yaml"):
        nuisance_parameters = {}

        super(GammapyLike, self).__init__(name, nuisance_parameters=nuisance_parameters)
        instrument = "veritas"
        if instrument == "veritas":
             config_file = yaml.load(open(config_file), Loader=yaml.FullLoader)
             self.data_dir = config_file['data']['anasum']
             self.output_dir = config_file['fileio']['outdir']
             self.on_region_radius = Angle("{} deg".format(np.sqrt(config_file['cuts']['th2cut'])))
             self.emin = config_file['selection']['emin']
             self.emax = config_file['selection']['emax']
             self.nbin = config_file['selection']['nbin']
             self.exclusion_on = config_file['selection']['exc_on_region_radius']
             self.exc_radius = config_file['selection']['exc_radius']
             datastore = DataStore.from_dir(self.data_dir)
             self.obs_table = datastore.obs_table
             obs_ids = self.obs_table['OBS_ID']
             available_irf = ["aeff", "edisp"]
             observations = datastore.get_observations(obs_ids, required_irf=available_irf)
             RA=self.obs_table['RA_OBJ']
             DEC=self.obs_table['DEC_OBJ']
             RA_OBJ = RA[0]
             DEC_OBJ = DEC[0]
             target_position = SkyCoord(ra=RA_OBJ, dec=DEC_OBJ, unit="deg", frame="icrs")
             on_region_radius = Angle("{} deg".format(np.sqrt(0.008)))
             on_region = CircleSkyRegion(center=target_position, radius=on_region_radius)

             exclusion_mask = []

             reg0 = CircleSkyRegion(center=SkyCoord(RA_OBJ, DEC_OBJ, unit="deg", frame="icrs"), radius=self.exclusion_on * u.deg, )
             reg1 = CircleSkyRegion(center=SkyCoord(81.9087, 21.937, unit="deg", frame="icrs"), radius=self.exc_radius * u.deg, )
             reg2 = CircleSkyRegion(center=SkyCoord(82.6806, 22.4623, unit="deg", frame="icrs"), radius=self.exc_radius * u.deg, )
             reg3 = CircleSkyRegion(center=SkyCoord(83.4118, 20.4742, unit="deg", frame="icrs"), radius=self.exc_radius * u.deg, )
             reg4 = CircleSkyRegion(center=SkyCoord(84.1099, 21.9931, unit="deg", frame="icrs"), radius=self.exc_radius * u.deg, )
             reg5 = CircleSkyRegion(center=SkyCoord(84.4112, 21.1425, unit="deg", frame="icrs"), radius=self.exc_radius * u.deg, )
             reg6 = CircleSkyRegion(center=SkyCoord(84.8629, 21.7629, unit="deg", frame="icrs"), radius=self.exc_radius * u.deg, )
             reg7 = CircleSkyRegion(center=SkyCoord(85.4782, 23.3262, unit="deg", frame="icrs"), radius=self.exc_radius * u.deg, )
             reg8 = CircleSkyRegion(center=SkyCoord(85.5166, 22.6603, unit="deg", frame="icrs"), radius=self.exc_radius * u.deg, )

             skydir = target_position.galactic
             geom = WcsGeom.create(npix=(1000, 1000), binsz=0.005, skydir=skydir, proj="TAN", frame="icrs")
             exclusion_mask = ~geom.region_mask([reg0, reg1, reg2, reg3, reg4, reg5, reg6, reg7, reg8])

             energy_ax = MapAxis.from_energy_bounds(self.emin, self.emax, nbin=self.nbin,  unit="TeV", name="energy")
             energy_ax_true = MapAxis.from_energy_bounds(self.emin, self.emax, nbin=self.nbin, unit="TeV", name="energy_true")
             
             geom = RegionGeom.create(region=on_region, axes=[energy_ax])
             dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_ax_true)

             dataset_maker = SpectrumDatasetMaker(containment_correction=False, selection=["counts", "exposure", "edisp"])
             bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)

             safe_mask_masker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=10)


             datasets = Datasets()
             for obs_id, observation in zip(obs_ids, observations):
                 dataset = dataset_maker.run(
                     dataset_empty.copy(name=str(obs_id)), observation
                 )
                 dataset_on_off = bkg_maker.run(dataset, observation)
                 dataset_on_off = safe_mask_masker.run(dataset_on_off, observation)
                 datasets.append(dataset_on_off)

                 
             self.dataset_stacked = Datasets(datasets).stack_reduce(name="stacked")

             

    def _get_gammapy_instance(self, likelihood_model):
#        self.model = likelihood_model
        for point_source in list(likelihood_model.point_sources.values()):
            pivot_energy = point_source.spectrum.main.Powerlaw.parameters['piv'].value
            pivot_eunit = point_source.spectrum.main.Powerlaw.parameters['piv'].unit
            index = point_source.spectrum.main.Powerlaw.parameters['index'].value
            k_value = point_source.spectrum.main.Powerlaw.parameters['K'].value
            k_unit = point_source.spectrum.main.Powerlaw.parameters['K'].unit
            self.spectral_model = PowerLawSpectralModel(
                index=index, 
                amplitude=k_value * u.Unit(k_unit), 
                reference=1 * u.Unit(pivot_eunit)
            )

        


        model = SkyModel(spectral_model=self.spectral_model, name="{}".format(self.obs_table['OBJECT'][0]))
        return model
        #self.dataset_stacked.models = [self.model]
        #return self.dataset_stacked






    def set_model(self, likelihood_model_instance):
        """
        Set the model to be used in the joint minimization. Must be a LikelihoodModel instance.
        """

        # This will take a long time if it's the first time we run, as it will select the data,
        # produce livetime cube, expomap, source maps and so on

        self._likelihood_model = likelihood_model_instance

#        self._dataset_stacked = self._get_gammapy_instance(likelihood_model_instance)
        #self._update_model_in_fermipy( update_dictionary = True, force_update = True)

        
    def get_log_like(self):
        """
        Return the value of the log-likelihood with the current values for the
        parameters stored in the ModelManager instance
        """

        # Update all sources on the fermipy side
        #self.set_model(likelihood_model_instance)
        #self._update_model_in_fermipy() #   # I still dont understnad this but should be easy if you dig a little deeper
        #should be something like this
        #self._update_model_in_gammapy()

        # Get value of the log likelihood
        model = self._get_gammapy_instance(self._likelihood_model)
        self.dataset_stacked.models = [model]                                                                                                                                                              
        try:

            value = self.dataset_stacked.stat_sum()
            print(value)

        except:

            raise

#        return value #- logfactorial(int(self._gta.like.total_nobs()))
        return value - logfactorial(int(271))

    def inner_fit(self):
        """
        This is used for the profile likelihood. Keeping fixed all parameters in the
        LikelihoodModel, this method minimize the logLike over the remaining nuisance
        parameters, i.e., the parameters belonging only to the model for this
        particular detector. If there are no nuisance parameters, simply return the
        logLike value.
        """
        return self.get_log_like()

    def get_number_of_data_points(self):
        """
        This returns the number of data points that are used to evaluate the likelihood.
        For binned measurements, this is the number of active bins used in the fit. For
        unbinned measurements, this would be the number of photons/particles that are
        evaluated on the likelihood
        """

        return int(6)  
