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
        Set the model to be used in the joint minimization.
        Must be a Astromodels Model instance.
        """

        self._likelihood_model: Model = likelihood_model_instance
        self._likelihood_model_converted = AstromodelConverter(self._likelihood_model)

    def get_log_like(self) -> float:
        """
        Return the value of the log-likelihood with the current values for the
        parameters stored in the model instance
        """
        self._likelihood_model_converted._update_parameters(self._likelihood_model)
        self._stacked.models = self.gammapy_model

        return -self._stacked.stat_sum()

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
        return np.prod(self._stacked.counts.data.shape)

    @property
    def model(self):
        return self._likelihood_model

    @property
    def gammapy_model(self):
        return self._likelihood_model_converted.gammapy_models
