import astropy.units as u
import numpy as np
from astromodels.core.parameter import Parameter
from gammapy.modeling.models import PowerLawSpectralModel, SkyModel

from gammapy_plugin.utils.gammapy_parser import (
    parameter_to_gammapy_dict,
    parse_gammapy_model,
)


def test_parameter_to_gammapy_dict():
    compare_dict = {"value": 2, "unit": u.TeV, "min": 1.0, "max": 3.0, "frozen": False}
    para = Parameter(
        name="test_paramter",
        value=2.0,
        unit=u.TeV,
        min_value=1.0,
        max_value=3.0,
        free=True,
    )
    para_dict = parameter_to_gammapy_dict(para)
    for k in compare_dict.keys():
        assert para_dict[k] == compare_dict[k]


def test_parse_gammapy_model():
    pl = PowerLawSpectralModel(
        amplitude=2.3e-11 * u.Unit("TeV-1 cm-2 s-1"),
        index=2.3,
        reference=1 * u.TeV,
    )
    sm = SkyModel(name="crab", spectral_model=pl)
    astromodel_parameters = parse_gammapy_model(sm)
    for k, v in astromodel_parameters.items():
        pn = k.split(".")[-1]
        if k == "empty.crab.amplitude":
            assert v.is_normalization is True
        elif k in ["empty.crab.index", "empty.crab.reference"]:
            pass
        else:
            raise AttributeError(f"Unkown parameter {k}")
        assert v.value == pl.parameters[pn].value
        assert v.unit == pl.parameters[pn].unit
        if v.min_value is not None:
            assert v.min_value == pl.parameters[pn].min
        else:
            assert pl.parameters[pn].min is np.nan
        if v.max_value is not None:
            assert v.max_value == pl.parameters[pn].max
        else:
            assert pl.parameters[pn].max is np.nan
        assert v.free is not pl.parameters[pn].frozen
