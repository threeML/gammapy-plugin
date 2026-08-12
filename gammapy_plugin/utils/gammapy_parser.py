import logging

import astropy.units as u
import numpy as np
from astromodels.core.parameter import Parameter
from gammapy.modeling.models.core import ModelBase

__all__ = ["parameter_to_gammapy_dict", "parse_gammapy_model"]

log = logging.getLogger(__name__)

_conversion_list = {
    u.Unit("keV-1 cm-2 s-1"): [1e9, u.Unit("TeV-1 cm-2 s-1")],
    u.Unit("keV"): [1e9, u.Unit("TeV")],
}


def parameter_to_gammapy_dict(para: Parameter) -> dict:
    """Converts a astromodel parameter to a dict able to be read in as a
    gammapy parameter.

    :param para: astromodel paramter
    :type para: Parameter
    :return: dict
    """
    para_dict = {}

    if para.unit in _conversion_list.keys():
        factor = _conversion_list[para.unit][0]
        para_unit = _conversion_list[para.unit][1]
    else:
        factor = 1.0
        para_unit = para.unit
    para_dict["value"] = para.value * factor
    para_dict["unit"] = para_unit
    val = np.nan
    if para.min_value is not None:
        val = para.min_value * factor
    para_dict["min"] = val
    val = np.nan
    if para.max_value is not None:
        val = para.max_value * factor
    para_dict["max"] = val
    para_dict["frozen"] = not para.free
    para_dict["prior"] = ""
    return para_dict


def parse_gammapy_model(gp_model: ModelBase, dataset_name: str = "empty") -> dict:
    """Returns dict of astromodels parameters with all parameters from the
    passed gammapy model.

    :param gp_model: gammapy model
    :param dataset_name: name of the dataset the model is associated to
        defaults to empty
    :type gp_model: ModelBase
    :type dataset_name: str
    """
    tp = gp_model.parameters.to_dict()
    parameters = {}
    for i in range(len(tp)):
        ttp = tp[i]
        norm = False
        if ttp["name"] in ["norm", "amplitude"]:
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
        log.debug("Created parameter")
        log.debug(
            parameters[dataset_name + "." + str(gp_model.name) + "." + ttp["name"]]
        )

    return parameters
