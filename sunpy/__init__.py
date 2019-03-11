"""
SunPy
=====

An open-source Python library for Solar Physics data analysis.

Web Links
---------
Homepage: http://sunpy.org
Documentation: http://docs.sunpy.org/en/stable/
"""
# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import os
import sys

__minimum_python_version__ = "3.6"


class UnsupportedPythonError(Exception):
    pass


if sys.version_info < tuple((int(val) for val in __minimum_python_version__.split('.'))):
    raise UnsupportedPythonError("sunpy does not support Python < {}".format(__minimum_python_version__))

# this indicates whether or not we are in the package's setup.py
try:
    _SUNPY_SETUP_
except NameError:
    import builtins
    builtins._SUNPY_SETUP_ = False

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

if not _SUNPY_SETUP_:
    from sunpy.util.config import load_config, print_config
    from sunpy.util import system_info
    from sunpy.tests.runner import SunPyTestRunner

    self_test = SunPyTestRunner.make_test_runner_in(os.path.dirname(__file__))

    # Load user configuration
    config = load_config()

    __all__ = ['config', 'self_test', 'system_info']
