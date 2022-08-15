__author__ = "J Michael Burgess"
import os
import uuid
from typing import List

import astropy.units as u
import gammapy.modeling.models.spectral as gpyspec
from astromodels import Function1D, Model


def _generate_uuid() -> str:
    """
    Generate a unique identifier for this function.

    :return: the UUID
    """
    return str(uuid.UUID(bytes=os.urandom(16), version=4)).replace("-", "_")


class SpectralModelGenerator:
    def __init__(self, function: Function1D) -> None:

        self._function: Function1D = function
        self._build_class()

    @property
    def class_def(self) -> gpyspec.SpectralModel:
        return self._class_def

    @property
    def class_name(self) -> str:
        return self._class_name

    def _build_class(self) -> None:

        self._class_name: str = f"Model{_generate_uuid()}"

        def evaluate(self, *args, **kwargs):

            return self._astro_func(args[0])

        class_dict = {}

        for i, (k, v) in enumerate(self._function.parameters.items()):

            is_norm = False

            if i == 0:
                is_norm = True

            class_dict[k] = gpyspec.Parameter(
                k, v.value * v._unit, is_norm=is_norm
            )

        class_dict["_astro_func"] = self._function
        class_dict["evaluate"] = evaluate

        self._class_def = type(
            self._class_name, (gpyspec.SpectralModel,), class_dict
        )


class GammapyModelWrapper:
    def __init__(self, model: Model) -> None:

        self._model: Model = model

        self._point_sources: List[gpyspec.SpectralModel] = []

        self._extended_sources = None

        # check for point sources

        self._build_point_sources()

    @property
    def point_sources(self) -> List[gpyspec.SpectralModel]:
        return self._point_sources

    def _build_point_sources(self) -> None:

        for i, (name, ps) in enumerate(self._model.point_sources.items()):

            # generate the class

            smc = SpectralModelGenerator(ps.spectrum.main.shape)

            self._point_sources.append(smc.class_def())
