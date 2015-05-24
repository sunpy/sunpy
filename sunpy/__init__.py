"""
SunPy
=====

An open-source Python library for Solar Physics data analysis.

Web Links
---------
Homepage: http://www.sunpy.org
Documentation: http://docs.sunpy.org
"""
from __future__ import absolute_import

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

try:
    from .version import githash as __githash__
except ImportError:
    __githash__ = ''

from sunpy.util.config import load_config, print_config
from sunpy.util import system_info
from sunpy.tests import main as self_test

# Load user configuration
config = load_config()