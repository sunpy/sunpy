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

    import logging

    # Use the root logger as a dummy log before initializing Astropy's logger
    log = logging.getLogger()

    from astropy.logger import _init_log

    #  read and translate settings for Astropy's logger
    if config.has_section('logger'):
        from astropy.logger import Conf as LoggerConf
        conf = LoggerConf()
        option_list = ['log_level', 'use_color', 'log_warnings', 'log_exceptions', 'log_to_file', 'log_file_path',
                       'log_file_level', 'log_file_format']
        for this_option in option_list:
            if config.has_option('logger', this_option):
                setattr(conf, this_option, config.get('logger', this_option))

    log = _init_log()

    __all__ = ['config', 'self_test', 'system_info']
