import numpy as np
from astromodels.core.parameter import Parameter
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
