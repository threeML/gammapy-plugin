from typing import Optional
import numpy as np
from astromodels import Model
from threeML.io.logging import setup_logger
from threeML.plugin_prototype import PluginPrototype
from gammapy_plugin.gammapy_converter import AstromodelConverter
from gammapy.datasets import Datasets, Dataset

log = setup_logger(__name__)

__instrument_name = "Gammapy"


class GammapyLike(PluginPrototype):
    def __new__(cls, *args, **kwargs):
        instance = object.__new__(cls)
        return instance

    def __init__(self, name: str, config_file: str = "config.yaml") -> None:
        """ """
        nuisance_parameters = {}

        super(GammapyLike, self).__init__(name, nuisance_parameters=nuisance_parameters)

    def set_datasets(self, datasets, **kwargs):
        """
        Set the Gammapy Dataset
        :param datasets: list of Gammapy datasets or a single Datasets object
        """
        if isinstance(datasets, list):
            self._datasets = Datasets()
            for d in datasets:
                self._datasets.append(d)
        elif isinstance(datasets, Datasets):
            self._datasets = datasets
        else:
            msg = "Datasets has to be list of Dataset or a single Datasets object"
            raise TypeError(msg)
        self._stacked = self._datasets.stack_reduce(name="stacked")

    def set_model(self, likelihood_model_instance: Model) -> None:
        """
        Set the model to be used in the joint minimization. Must be a LikelihoodModel instance.
        """

        # This will take a long time if it's the first time we run, as it will select the data,
        # produce livetime cube, expomap, source maps and so on

        self._likelihood_model: Model = likelihood_model_instance

        self._gammapy_wrapper = GammapyModelWrapper(self._likelihood_model)

        self._gammapy_model = SkyModel(
            spectral_model=self._gammapy_wrapper.point_sources[0],
            name=f"{self.obs_table['OBJECT'][0]}",
        )

    #        self._dataset_stacked = self._get_gammapy_instance(likelihood_model_instance)
    # self._update_model_in_fermipy( update_dictionary = True, force_update = True)

    @property
    def gammapy_model(self) -> Optional[SkyModel]:
        return self._gammapy_model

    @property
    def gammapy_wrapper(self) -> GammapyModelWrapper:
        return self._gammapy_wrapper

    def get_log_like(self) -> float:
        """
        Return the value of the log-likelihood with the current values for the
        parameters stored in the ModelManager instance
        """

        # Update all sources on the fermipy side
        # self.set_model(likelihood_model_instance)
        # self._update_model_in_fermipy() #   # I still dont understnad this but should be easy if you dig a little deeper
        # should be something like this
        # self._update_model_in_gammapy()

        # Get value of the log likelihood

        # model = self._get_gammapy_instance(self._likelihood_model)

        # maybe do this onley once?

        self.dataset_stacked.models = [self.gammapy_model]

        try:

            value = self.dataset_stacked.stat_sum()
        #            print(value)

        except:

            raise

        #        return value #- logfactorial(int(self._gta.like.total_nobs()))
        return -value  # - logfactorial(int(271))

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

        return (
            self.dataset_stacked.data_shape[0]
            * self.dataset_stacked.data_shape[1]
            * self.dataset_stacked.data_shape[2]
        )
