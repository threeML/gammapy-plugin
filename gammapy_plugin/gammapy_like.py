import logging
from pathlib import Path
from typing import TYPE_CHECKING

import gammapy
import numpy as np
from astromodels.core.model import Model
from astromodels.functions.priors import Truncated_gaussian
from gammapy.datasets import Dataset, Datasets
from gammapy.modeling.models import DatasetModels, ModelBase, Models
from packaging.version import Version
from threeML.plugin_prototype import PluginPrototype

from gammapy_plugin.converter import AstromodelConverter
from gammapy_plugin.io.plot_model import plot_model
from gammapy_plugin.utils.gammapy_parser import (
    parameter_to_gammapy_dict,
    parse_gammapy_model,
)

if TYPE_CHECKING:
    from threeML.io.plotting.data_residual_plot import ResidualPlot


__all__ = ["GammapyLike"]

log = logging.getLogger(__name__)

__instrument_name = "gammapy"

gammapy_version = Version(gammapy.__version__)


class GammapyLike(PluginPrototype):
    """A plugin for including Gammapy datasets."""

    def __new__(cls, *args, **kwargs) -> PluginPrototype:
        instance = object.__new__(cls)
        return instance

    def __init__(self, name: str, **kwargs) -> None:
        """
        Initialize the GammapyLike plugin.

        :param name: Name of the plugin
        :type name: str
        """
        nuisance_parameters = kwargs.get("nuisance_parameters", {})
        super(GammapyLike, self).__init__(name, nuisance_parameters=nuisance_parameters)
        self._frame: str = kwargs.get("frame", "icrs")
        self._sources: list[str] | None = kwargs.get("sources", None)
        self._nuisance_mapping = {}
        self._background_models: dict = kwargs.get("background_models", {})
        self._background_models_mapper: dict = {}
        self._nuisance_parameters_dicts: dict = {}
        if len(self._background_models.keys()) > 0:
            self._parse_background_models()

    def set_datasets(
        self,
        datasets: Dataset | Datasets | list[Dataset] | Path | str,
        mode: str = "individual",
        stacked_name: str = "stacked",
    ) -> None:
        """Set the Gammapy Dataset.

        :param datasets: list of Gammapy datasets or a single Dataset object
        :param mode: individual or stacked - defaults to individual, stacked
            stacks the passed datasets
        :param stacked_name: name of the stacked datasset if mode is stacked
        """
        if mode not in [
            "individual",
            "stacked",
        ]:
            raise ValueError("mode needs to be individual or stacked")
        if isinstance(datasets, (str, Path)):
            log.warn(
                "You have provided a path to a file - we will assume this is a Datasets"
                " file. If not this will likely fail"
            )
            datasets = Datasets.read(datasets)
            self._datasets = datasets
            if mode == "stacked":
                self._datasets = Datasets(
                    self._datasets.stack_reduce(name=stacked_name)
                )

        if isinstance(datasets, list):
            self._datasets = Datasets()
            for d in datasets:
                self._datasets.append(d)
            if mode == "stacked":
                self._datasets = Datasets(
                    self._datasets.stack_reduce(name=stacked_name)
                )

        elif isinstance(datasets, Datasets):
            self._datasets = datasets
            if mode == "stacked":
                self._datasets = Datasets(
                    self._datasets.stack_reduce(name=stacked_name)
                )
        elif isinstance(datasets, Dataset):
            self._datasets = Datasets(datasets)
            if mode == "stacked":
                log.info("Only using a single dataset - can not stack that")
        else:
            msg = "datasets has to be list of Dataset,"
            msg += " a single Datasets or Dataset object"
            raise TypeError(msg)

    def set_sources(self, sources: list[str] | str = None) -> None:
        """
        Set the sources to be used by this plugin - No need to specify bkg models

        :param sources: Source(s) to be used in the analysis defaults to all
        :type sources: list of str or str
        """

        if isinstance(sources, list):
            self._sources = sources
        elif isinstance(sources, str):
            self._sources = [sources]
        elif sources is None:
            self._sources = None
        else:
            raise ValueError("")

    def set_model(
        self,
        likelihood_model: Model,
        converted_model: AstromodelConverter = None,
    ) -> None:
        """
        Set the model to be used in the joint minimization.

        :param likelihood_model: astromodels model
        :param converted_model: converted astromodels
        :type likelihood_model: Model
        :type converted_model: AstromodelConverter
        """

        if self._sources is None:
            log.debug(
                "If you want to specify sources for this Plugin you MUST do so before"
            )
        else:
            log.debug(f"Will use {self._sources} for this plugin")

        self._likelihood_model: Model = likelihood_model
        if converted_model is not None:
            self._likelihood_model_converted: AstromodelConverter = converted_model
        elif hasattr(self, "_likelihood_model_converted"):
            pass
        else:
            self._likelihood_model_converted: AstromodelConverter = AstromodelConverter(
                model=self._likelihood_model, frame=self._frame
            )

        self._update_gammapy_model_list()
        self._assign_models()
        self._set_get_log_like()

    def _update_gammapy_model_list(self) -> Models:
        """
        Update the list of gammapy models.
        """

        # this will only be run after setting the model
        if hasattr(self, "_likelihood_model_converted"):
            if self._sources is not None:
                tmp2 = [
                    x
                    for x in self._likelihood_model_converted.gammapy_models
                    if x.name in self._sources
                ]

                tmp = [*tmp2]
            else:
                tmp2 = [x for x in self._likelihood_model_converted.gammapy_models]
                tmp = [*tmp2]
        else:
            tmp = []
            tmp2 = []
        self._global_models = Models(tmp2)
        if hasattr(self, "_background_models"):
            for m in list(self._background_models.values()):
                tmp.append(m)
        self._gammapy_model = Models(tmp)

    def _assign_models(self):
        """
        Assign the gammapy models to the datasets
        """
        for d in self._datasets:
            d.models = []
            tmp = []
            for g in self._global_models:
                tmp.append(g)
            if d.name in self._background_models_mapper.keys():
                tmp.append(
                    self._background_models[self._background_models_mapper[d.name]]
                )
            d.models = Models(tmp)

    def set_background_models(
        self, bkg_model: ModelBase | list | Models | DatasetModels
    ) -> None:
        """
        Set the gammapy background models (e.g. FoVBackgroundModel)

        :param bkg_model: Background model(s) :type bkg_model: ModelBase or list of
            ModelBase or Models or DatasetModels.
        """
        if isinstance(bkg_model, ModelBase):
            bkg_model = [bkg_model]
        else:
            if not isinstance(bkg_model, (list, Models, DatasetModels)):
                raise TypeError(
                    "either pass a singular gammapy model or a list of models"
                )
        for b in bkg_model:
            self._background_models[b.name] = b
            self._background_models_mapper[b.datasets_names[0]] = b.name
        self._parse_background_models()
        self._update_gammapy_model_list()
        self._assign_models()

    def _parse_background_models(self):
        """
        Parse the background models and link the gammapy parameters to
        nuissance parameters of this plugin and set the prior.
        """
        # TODO: way of manually specifying the priors
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
        """
        Update the background models with the current values of the nuisance
        parameters.
        """
        # TODO: rewrite this to directly update
        for k, v in self._nuisance_parameters.items():
            self._nuisance_parameters_dicts[k]["value"] = v.value
            p = self._nuisance_mapping[k]
            # TODO find an elegant way to not hardcode this
            self._background_models[p[1]].parameters[p[2]].update_from_dict(
                self._nuisance_parameters_dicts[k]
            )

    def _set_get_log_like(self):
        if hasattr(self._datasets, "_stat_sum_likelihood"):
            func = self._datasets._stat_sum_likelihood
        elif hasattr(self._datasets, "stat_sum_likelihood"):
            func = self._datasets.stat_sum_likelihood
        else:
            raise AttributeError(
                "gammapy.Datasets has neither _stat_sum_likelihood nor "
                "stat_sum_likelihood - something went fundamentally wrong!"
            )

        self.log_like_func = func

    def get_log_like(self) -> float:
        """
        Return the value of the log-likelihood with the current values for
        the parameters stored in the model instance.
        """
        self._likelihood_model_converted._update_parameters()
        self._update_background_models()
        return -0.5 * self.log_like_func()

    def inner_fit(self):
        return self.get_log_like()

    def get_number_of_data_points(self) -> np.int64:
        """
        Return the number of data points in the datasets.
        """
        return np.sum([np.prod(d.counts.data.shape) for d in self._datasets])

    @property
    def datasets(self) -> Datasets:
        """Gammapy datasets of the plugin."""
        return self._datasets

    @property
    def model(self) -> Model:
        """Astromodels model of the plugin."""
        return self._likelihood_model

    @property
    def astromodel_converter(self) -> AstromodelConverter:
        """AstromodelConverter object used for this plugin."""
        return self._likelihood_model_converted

    @property
    def gammapy_model(self) -> Models:
        """List of all the Gammapy SkyModels."""

        if not hasattr(self, "_gammapy_model"):
            self._update_gammapy_model_list()
        return self._gammapy_model

    @property
    def frame(self) -> str:
        """Coordinate Frame of the plugin."""
        return self._frame

    def display_model(
        self,
        *args,
        **kwargs,
    ) -> "ResidualPlot":

        return plot_model(self._datasets, *args, **kwargs)
