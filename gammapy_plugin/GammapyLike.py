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
    def __new__(cls, *args, **kwargs) -> PluginPrototype:
        instance = object.__new__(cls)
        return instance

    def __init__(self, name: str, **kwargs) -> None:
        nuisance_parameters = {}
        super(GammapyLike, self).__init__(name, nuisance_parameters=nuisance_parameters)
        self._frame = kwargs.get("frame", "icrs")
        self._sources = kwargs.get("sources", None)

    def set_datasets(self, datasets, **kwargs):
        """
        Set the Gammapy Dataset
        :param datasets: list of Gammapy datasets or a single Datasets object
        """
        if isinstance(datasets, list):
            self._datasets = Datasets()
            for d in datasets:
                self._datasets.append(d)
            self._stacked = self._datasets.stack_reduce(name="stacked")
        elif isinstance(datasets, Datasets):
            self._datasets = datasets
            self._stacked = self._datasets.stack_reduce(name="stacked")
        elif isinstance(datasets, Dataset):
            self._datasets = Datasets(datasets)
            self._stacked = datasets.copy(name="stacked")
        else:
            msg = "Datasets has to be list of Dataset or a single Datasets object"
            raise TypeError(msg)

    def set_sources(self, sources: list = None):
        """
        Set the sources to be used by this plugin
        Needed e.g. for assigning different background models to the plugins
        """
        # todo assert source is in model
        assert isinstance(sources, list) or isinstance(
            sources, None
        ), "Wrong source type"
        self._sources = sources

    def set_model(
        self, likelihood_model_instance: Model, gammapy_model: ModelBase = None
    ) -> None:
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
        self._stacked.models = self.gammapy_model

        return -self._stacked._stat_sum_likelihood()

    def inner_fit(self):
        return self.get_log_like()

    def get_number_of_data_points(self):
        return np.prod(self._stacked.counts.data.shape)

    @property
    def stacked(self) -> Dataset:
        return self._stacked

    @property
    def model(self) -> Model:
        return self._likelihood_model

    @property
    def gammapy_model(self) -> list[SkyModel]:
        return self._likelihood_model_converted.gammapy_models

    @property
    def frame(self) -> str:
        return self._frame
