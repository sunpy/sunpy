"""
SunPy
=====

An open-source Python library for Solar Physics data analysis.

Web Links
---------
Homepage: http://www.sunpy.org
Documentation: http://sunpy.readthedocs.org/en/latest/index.html
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

# Sample data
from sunpy.data.sample import (AIA_171_IMAGE, RHESSI_IMAGE, EIT_195_IMAGE,
                               RHESSI_EVENT_LIST, CALLISTO_IMAGE)

# Load user configuration
config = load_config()
