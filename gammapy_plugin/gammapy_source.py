import collections
from typing import Optional
import numpy as np
from astromodels.sources.source import Source, SourceType
from astromodels.core.tree import Node
from astromodels.core.parameter import Parameter
from astromodels.utils.pretty_list import dict_to_list
from astromodels.functions.function import Function, FunctionMeta
from astromodels.functions.functions_1D import Constant
from gammapy.modeling.models.core import ModelBase


def parameter_to_gammapy_dict(para: Parameter) -> dict:
    para_dict = {}
    para_dict["value"] = para.value
    para_dict["unit"] = para.unit
    val = np.nan
    if para.min_value is not None:
        val = para.min_value
    para_dict["min"] = val
    val = np.nan
    if para.max_value is not None:
        val = para.max_value
    para_dict["max"] = val
    para_dict["frozen"] = not para.free
    para_dict["prior"] = ""
    return para_dict


def parse_gammapy_model(gp_model: ModelBase, dataset_name: str) -> dict:
    """
    Returns astromodels dict with all parameters from the
    passed gammapy model
    """
    tp = gp_model.parameters.to_dict()
    parameters = {}
    for i in range(len(tp)):
        ttp = tp[i]
        norm = False
        if ttp["name"] == "norm":
            # TODO check if this covers all cases
            norm = True
        min_v = ttp["min"]
        max_v = ttp["max"]
        if np.isnan(min_v):
            min_v = None
        if np.isnan(max_v):
            max_v = None
        parameters[dataset_name + "." + str(gp_model.name) + "." + ttp["name"]] = (
            Parameter(
                name=dataset_name + "_" + str(gp_model.name) + "_" + ttp["name"],
                value=ttp["value"],
                min_value=min_v,
                max_value=max_v,
                is_normalization=norm,
                free=not ttp["frozen"],
                unit=ttp["unit"],
                desc=f"Gammapy Model Parameter {ttp['name']}",
            )
        )
    return parameters


class GammapySource(Source, Node):
    # TODO does not work ... issue with astromodels not knowing the
    # used funtcion
    """
    A dummy source for using a GammaPy Model

    :param name: name of the source
    :param gammapy_model: the Gammapy Model
    """

    def __init__(self, name: str, gammapy_model: ModelBase) -> None:
        self._gammapy_model = gammapy_model
        Source.__init__(self, [], SourceType.EXTENDED_SOURCE)
        Node.__init__(self, name)
        self._get_parameters()
        fct = {"description": "GP Model"}
        fct = self._add_paras_fct_description(fct)
        gp_shape = GPFunction(
            name=name, function_definition=fct, parameters=self._parameters
        )

        # needed for all astromodel display functionalities to work
        # spectrum_node = Node("spectrum")
        # self._add_child(spectrum_node)
        self._add_child(gp_shape)

    def _add_paras_fct_description(self, fct: dict) -> dict:
        fct["parameters"] = {}
        for k, v in self.parameters.items():
            fct["parameters"][k] = v.to_dict()
        return fct

    def _get_parameters(self) -> None:
        """
        Load all the parameters from the GP model and create them in the
        astromodel source
        """
        tp = self._gammapy_model.parameters.to_dict()
        parameters = {}
        for i in range(len(tp)):
            ttp = tp[i]
            norm = False
            if ttp["name"] == "norm":
                # TODO check if this covers all cases
                norm = True
            min_v = ttp["min"]
            max_v = ttp["max"]
            if np.isnan(min_v):
                min_v = None
            if np.isnan(max_v):
                max_v = None
            parameters[ttp["name"]] = Parameter(
                name=ttp["name"],
                value=ttp["value"],
                min_value=min_v,
                max_value=max_v,
                is_normalization=norm,
                free=not ttp["frozen"],
                unit=ttp["unit"],
                desc=f"Gammapy Model Parameter {ttp['name']}",
            )
        self._parameters = parameters

    @property
    def parameters(self) -> dict:
        return self._parameters

    @property
    def free_parameters(self) -> list[Parameter]:
        return [p for k, p in self._parameters.items() if p.free]

    @property
    def has_free_parameters(self) -> bool:
        ret = False
        for k, v in self._parameters.items():
            if v.free:
                ret = True
                break
        return ret

    @property
    def gammapy_model(self) -> ModelBase:
        """
        Gammapy model set when initiating
        """
        return self._gammapy_model

    def _repr__base(self, rich_output: bool = False) -> list:
        # TODO
        repr_dict = collections.OrderedDict()
        key = "%s (gammapy source)" % self.name
        repr_dict[key] = collections.OrderedDict()
        for component_name, component in list(self.components.items()):
            repr_dict[key][component_name] = component.to_dict(minimal=True)

        return dict_to_list(repr_dict, rich_output)


class GPFunction(Function):
    r"""
    description :

    parameters :

        norm :

            desc : Integral between a and b
            initial value : 1
            is_normalization : True
            transformation : log10
            min : 1e-30
            max : 1e3
            delta : 0.1

    """

    def evaluate(self, x, norm):
        pass

    def _set_units(self, x_unit, y_unit):
        pass
