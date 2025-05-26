import numpy as np
from astromodels.core.parameter import Parameter
from gammapy.modeling.models.core import ModelBase
from threeML.io.logging import setup_logger

log = setup_logger(__name__)


def parameter_to_gammapy_dict(para: Parameter) -> dict:
    """
    Converts a astromodel parameter to a dict able to be read
    in as a gammapy parameter

    :param para: astromodel paramter
    :type para: Parameter

    :return: dict
    """
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


def parse_gammapy_model(gp_model: ModelBase, dataset_name: str = "empty") -> dict:
    """
    Returns dict of astromodels parameters with all parameters from the
    passed gammapy model

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
        log.debug("Created parameter")
        log.debug(
            parameters[dataset_name + "." + str(gp_model.name) + "." + ttp["name"]]
        )

    return parameters
