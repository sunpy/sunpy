"""
SunPy
=====

An open-source Python library for Solar Physics data analysis.

Web Links
---------
Homepage: http://sunpy.org
Documentation: http://docs.sunpy.org/en/stable/
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

import os
from sunpy.util.config import load_config, print_config
from sunpy.util import system_info
from sunpy.tests.runner import SunPyTestRunner

self_test = SunPyTestRunner.make_test_runner_in(os.path.dirname(__file__))

# Load user configuration
config = load_config()

__all__ = ['config', 'self_test', 'system_info']
