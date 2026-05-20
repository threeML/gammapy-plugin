"""Top-level package for Gammapy Plugin."""

from ._version import get_versions

from .gammapy_like import GammapyLike

__version__ = get_versions()["version"]
del get_versions


__all__ = ["GammapyLike"]
