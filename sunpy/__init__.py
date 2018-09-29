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

try:
    _ASTROPY_SETUP_
except NameError:
    _ASTROPY_SETUP_ = False


if not _ASTROPY_SETUP_:
    from sunpy.util.config import load_config, print_config
    from sunpy.util import system_info
    from sunpy.tests import main as self_test

    
    # Load user configuration
    config = load_config()

    __all__ = ['config', 'self_test', 'system_info']
