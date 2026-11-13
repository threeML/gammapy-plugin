import copy

import pytest
from astromodels.core.model import Model
from astromodels.functions import Powerlaw
from astromodels.sources import PointSource

from gammapy_plugin.converter import AstromodelConverter


def test_gammapy_updates_original_astromodels():
    pl = Powerlaw()
    ps = PointSource(
        "test",
        ra=0,
        dec=0,
        spectral_shape=pl,
    )

    am_model = Model(ps)
    converter = AstromodelConverter(am_model)

    gp_model = converter.gammapy_models[0]

    name = pl.index.path
    gp_parameter = gp_model.parameters[name]

    # Check that the actual parameter is linked.
    assert gp_parameter._am_parameter is pl.index

    # Direct Gammapy value update.
    gp_parameter.value = -2.5

    assert pl.index.value == pytest.approx(-2.5)
    assert am_model.parameters[name].value == pytest.approx(-2.5)

    # Optimizer-style update via factor.
    gp_parameter.factor = gp_parameter.transform(-2.8)

    assert pl.index.value == pytest.approx(-2.8)


def test_gammapy_model_copy():
    pl = Powerlaw()
    ps = PointSource(
        "test",
        ra=0,
        dec=0,
        spectral_shape=pl,
    )

    am_model = Model(ps)
    converter = AstromodelConverter(am_model)

    gp_model = converter.gammapy_models[0]
    copied_model = gp_model.copy()

    name = pl.index.path

    gp_parameter = copied_model.parameters[name]
    copied_function = copied_model.spectral_model._astromodel_function
    copied_am_parameter = copied_function.parameters["index"]

    # Check the binding itself!
    assert gp_parameter._am_parameter is copied_am_parameter

    # The copy must be independent of the original.
    assert copied_am_parameter is not pl.index

    # Changing the copied parameter updates its
    # own astromodels function.
    gp_parameter.value = -3.0

    assert copied_am_parameter.value == pytest.approx(-3.0)

    # The original astromodels model stays unchanged.
    assert pl.index.value == pytest.approx(-2.01)

    # Also check the original Gammapy model.
    assert gp_model.parameters[name].value == pytest.approx(-2.01)


@pytest.mark.parametrize(
    "copy_method",
    [
        copy.deepcopy,
        lambda model: model.copy(),
    ],
    ids=["python-deepcopy", "gammapy-copy"],
)
def test_copy_preserves_astromodels_bindings(copy_method):
    pl = Powerlaw()
    ps = PointSource("test", ra=0, dec=0, spectral_shape=pl)
    am_model = Model(ps)

    converter = AstromodelConverter(am_model)
    gp_model = converter.gammapy_models[0]

    original_value = pl.index.value
    copied_model = copy_method(gp_model)

    gp_par = copied_model.spectral_model.parameters[pl.index.path]
    am_par = copied_model.spectral_model._astromodel_function.parameters["index"]

    # Both parameters must be the same underlying object.
    assert gp_par._am_parameter is am_par

    # The copy must be independent.
    assert am_par is not pl.index

    # Setting the copied parameter propagates.
    gp_par.value = -3.0

    assert am_par.value == pytest.approx(-3.0)

    # The original model remains unchanged.
    assert pl.index.value == pytest.approx(original_value)
    assert gp_model.spectral_model.parameters[pl.index.path].value == pytest.approx(
        original_value
    )
