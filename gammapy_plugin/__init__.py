# -*- coding: utf-8 -*-
"""Top-level package for Gammapy Plugin."""

# please replace my name as this is boiler plate from
# my cookie cutter

# __author__ = """J. Michael Burgess"""
# __email__ = 'jburgess@mpe.mpg.de'


# from .utils.logging import gammapy_plugin_config, show_configuration

from . import _version

__version__ = _version.get_versions()["version"]
