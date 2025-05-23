from importlib.metadata import version
from typing import Union

import numpy as np
from astromodels.core.model import Model
from astromodels.functions.priors import Truncated_gaussian
from gammapy.datasets import Dataset, Datasets
from gammapy.modeling.models import DatasetModels, ModelBase, Models, SkyModel
from threeML.io.logging import setup_logger
from threeML.plugin_prototype import PluginPrototype

from gammapy_plugin.converter import AstromodelConverter
from gammapy_plugin.gammapy_source import (parameter_to_gammapy_dict,
                                           parse_gammapy_model)

GAMMAPY_VERSION = version("gammapy")
GAMMAPY_VERSION_MAJOR = int(GAMMAPY_VERSION.split(".")[0])

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
        nuisance_parameters = kwargs.get("nuisance_parameters", {})
        super(GammapyLike, self).__init__(name, nuisance_parameters=nuisance_parameters)
        self._frame = kwargs.get("frame", "icrs")
        self._sources = kwargs.get("sources", None)
        self._nuisance_mapping = {}
        self._background_models = kwargs.get("background_models", {})
        self._nuisance_parameters_dicts = {}
        if len(self._background_models.keys()) > 0:
            self._parse_background_models()

    def set_datasets(
        self,
        datasets: Union[Dataset, Datasets, list[Dataset]],
        mode: str = "individual",
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
            msg = "datasets has to be list of Dataset,"
            msg += " a single Datasets or Dataset object"
            raise TypeError(msg)

    def set_sources(self, sources: list = None) -> None:
        """
        Set the sources to be used by this plugin
        Needed e.g. for assigning different background models to individual plugins
        """
        # todo assert source is in model
        assert isinstance(sources, list) or sources is None, "Wrong source type"
        self._sources = sources

    def set_model(
        self,
        likelihood_model_instance: Model,
        converted_model: AstromodelConverter = None,
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
        if converted_model is not None:
            self._likelihood_model_converted = converted_model
        elif hasattr(self, "_likelihood_model_converted"):
            pass
        else:
            self._likelihood_model_converted = AstromodelConverter(
                model=self._likelihood_model, frame=self._frame
            )
        self._set_gammapy_model()
        for d in self._datasets:
            d.models = self.gammapy_model

    def set_background_models(self, bkg_model):
        if isinstance(bkg_model, ModelBase):
            bkg_model = [bkg_model]
        else:
            assert isinstance(
                bkg_model, (list, Models, DatasetModels)
            ), "either pass a singular gammapy model or a list of models"
        for b in bkg_model:
            self._background_models[b.name] = b
        self._parse_background_models()
        self._set_gammapy_model()
        for d in self._datasets:
            d.models = self.gammapy_model

    def _parse_background_models(self):
        """ """
        for name, bkg in self._background_models.items():
            bkg_paras = parse_gammapy_model(bkg, self._name)
            for k, v in bkg_paras.items():
                para_path = k.split(".")
                self._nuisance_mapping[k] = para_path
                self._nuisance_parameters[k] = v
                if v.is_normalization and v.free:
                    self._nuisance_parameters[k].prior = Truncated_gaussian(
                        mu=1.0, sigma=0.1, lower_bound=0.2, upper_bound=1.8
                    )
                self._nuisance_parameters_dicts[k] = parameter_to_gammapy_dict(v)

    def _update_background_models(self):
        # TODO rewrite this to directly update
        for k, v in self._nuisance_parameters.items():
            self._nuisance_parameters_dicts[k]["value"] = v.value
            p = self._nuisance_mapping[k]
            # TODO find an elegant way to not hardcode this
            self._background_models[p[1]].parameters[p[2]].update_from_dict(
                self._nuisance_parameters_dicts[k]
            )

    def get_log_like(self) -> float:
        """
        Return the value of the log-likelihood with the current values for the
        parameters stored in the model instance
        """
        self._likelihood_model_converted._update_parameters()
        self._update_background_models()
        # TODO: maybe too costly
        if GAMMAPY_VERSION_MAJOR > 1:
            return -self._datasets._stat_sum_likelihood()
        else:
            return -self._datasets.stat_sum()

    def inner_fit(self):
        return self.get_log_like()

    def get_number_of_data_points(self) -> np.int64:
        # TODO: check if this works for all allowed data types
        return np.sum([np.prod(d.counts.data.shape) for d in self._datasets])

    def distribute_covariance(self):
        raise NotImplementedError

    @property
    def datasets(self) -> Dataset:
        """
        Gammapy dataset of the plugin
        """
        return self._datasets

    @property
    def model(self) -> Model:
        """
        Astromodels model of the plugin
        """
        return self._likelihood_model

    def _set_gammapy_model(self) -> Models:
        """
        List of all the Gammapy SkyModels
        """
        # TODO do not do that in the property
        if hasattr(self, "_likelihood_model_converted"):
            if self._sources is not None:
                tmp = [
                    x
                    for x in self._likelihood_model_converted.gammapy_models
                    if x.name in self._sources
                ]
            else:
                tmp = [x for x in self._likelihood_model_converted.gammapy_models]
        else:
            tmp = []
        if hasattr(self, "_background_models"):
            for m in self._background_models.values():
                tmp.append(m)
        self._gammapy_model = Models(tmp)

    @property
    def gammapy_model(self):
        if not hasattr(self, "_gammapy_model"):
            self._set_gammapy_model()
        return self._gammapy_model

    @property
    def frame(self) -> str:
        """
        Coordinate Frame of the plugin
        """
        return self._frame


# TODO: setter method for frame --> change all the frames
# --> may also needs to be done for astromodels!
