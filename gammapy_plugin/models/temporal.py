import logging

from astromodels.functions.function import Function
from gammapy.modeling.models import (
    TemporalModel,
)

log = logging.getLogger(__name__)


class TemporalModelConverted(TemporalModel):
    def __init__(self, function: Function) -> None:
        raise NotImplementedError("Check how this is handled in gammapy")
