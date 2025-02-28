from typing import Union
import numpy as np
from astromodels import Model
from threeML.io.logging import setup_logger
from threeML.plugin_prototype import PluginPrototype
from gammapy_plugin.gammapy_converter import AstromodelConverter
from gammapy.datasets import Datasets, Dataset
from gammapy.modeling.models import SkyModel

log = setup_logger(__name__)

__instrument_name = "Gammapy"


class GammapyLike(PluginPrototype):
    """
    A plugin for including instruments supported by Gammapy
    After initiating you need to set_datasets() to add Gammapy dataset(s) to the
    Plugin

    """

    def __new__(cls, *args, **kwargs) -> PluginPrototype:
        instance = object.__new__(cls)
        return instance

    def __init__(self, name: str, **kwargs) -> None:
        nuisance_parameters = {}
        super(GammapyLike, self).__init__(name, nuisance_parameters=nuisance_parameters)
        self._frame = kwargs.get("frame", "icrs")
        self._sources = kwargs.get("sources", None)

    def set_datasets(
        self, datasets: Union[Dataset, Datasets], mode: str = "individual", **kwargs
    ) -> None:
        """
        Set the Gammapy Dataset
        :param datasets: list of Gammapy datasets or a single Dataset object
        :param mode: individual or stacked - defaults to individual
        """
        assert mode in [
            "individual",
            "stacked",
        ], "mode needs to be individual or stacked"
        if isinstance(datasets, list):
            self._datasets = Datasets()
            for d in datasets:
                self._datasets.append(d)
            if mode == "stacked":
                self._datasets = Datasets(self._datasets.stack_reduce(name="stacked"))

        elif isinstance(datasets, Datasets):
            self._datasets = datasets
            if mode == "stacked":
                self._datasets = Datasets(self._datasets.stack_reduce(name="stacked"))
        elif isinstance(datasets, Dataset):
            self._datasets = Datasets(datasets)
            if mode == "stacked":
                log.info("Only using a single dataset - can not stack that")
        else:
            msg = "datasets has to be list of Dataset, a single Datasets or Dataset object"
            raise TypeError(msg)

    def set_sources(self, sources: list = None) -> None:
        """
        Set the sources to be used by this plugin
        Needed e.g. for assigning different background models to individual plugins
        """
        # todo assert source is in model
        assert isinstance(sources, list) or sources is None, "Wrong source type"
        self._sources = sources

    def set_model(self, likelihood_model_instance: Model) -> None:
        """
        Set the model to be used in the joint minimization.
        Must be a Astromodels Model instance.
        """

        if self._sources is None:
            log.warning(
                "If you want to specify sources for this Plugin you MUST do so before"
            )
        else:
            log.info(f"Will use {self._sources} for this plugin")
        self._likelihood_model: Model = likelihood_model_instance
        self._likelihood_model_converted = AstromodelConverter(
            self._likelihood_model, self._frame, self._sources
        )

    def get_log_like(self) -> float:
        """
        Return the value of the log-likelihood with the current values for the
        parameters stored in the model instance
        """
        self._likelihood_model_converted._update_parameters()
        for d in self._datasets:
            d.models = [*self.gammapy_model]

        return -self._datasets._stat_sum_likelihood()

    def inner_fit(self):
        return self.get_log_like()

    def get_number_of_data_points(self) -> np.int64:
        return np.sum([np.prod(d.counts.data.shape) for d in self._datasets])

    @property
    def dataset(self) -> Dataset:
        """
        Gammapy dataset of the plugin
        """
        return self._dataset

    @property
    def model(self) -> Model:
        """
        Astromodels model of the plugin
        """
        return self._likelihood_model

    @property
    def gammapy_model(self) -> list[SkyModel]:
        """
        List of all the Gammapy SkyModels
        """
        return [*self._likelihood_model_converted.gammapy_models]

    @property
    def frame(self) -> str:
        """
        Coordinate Frame of the plugin
        """
        return self._frame
