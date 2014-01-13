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

import warnings

__version__ = '0.3.3'

try:
    from sunpy.version import version as __version__
except ImportError:
    warnings.warn('Missing version.py; you need to run setup.py', Warning)

from sunpy.util.config import load_config, print_config
from sunpy.util import system_info
from sunpy.tests import main as self_test

# Sample data
from sunpy.data.sample import (AIA_171_IMAGE, RHESSI_IMAGE, EIT_195_IMAGE,
                               RHESSI_EVENT_LIST, CALLISTO_IMAGE)

# Load user configuration
config = load_config()
