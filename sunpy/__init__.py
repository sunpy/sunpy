"""
SunPy
=====

An open-source Python library for Solar Physics data analysis.

Web Links
---------
Homepage: https://sunpy.org
Documentation: https://docs.sunpy.org/en/stable/
"""
# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import os
import sys

__minimum_python_version__ = "3.6"


class UnsupportedPythonError(Exception):
    pass


if sys.version_info < tuple(int(val) for val in __minimum_python_version__.split('.')):
    # This has to be .format to keep backwards compatibly.
    raise UnsupportedPythonError("Sunpy does not support Python < {}".format(__minimum_python_version__))

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

__citation__ = """\
@ARTICLE{2015CS&D....8a4009S,
   author = {{SunPy Community}, T. and {Mumford}, S.~J. and {Christe}, S. and
    {P{\'e}rez-Su{\'a}rez}, D. and {Ireland}, J. and {Shih}, A.~Y. and
    {Inglis}, A.~R. and {Liedtke}, S. and {Hewett}, R.~J. and {Mayer}, F. and
    {Hughitt}, K. and {Freij}, N. and {Meszaros}, T. and {Bennett}, S.~M. and
    {Malocha}, M. and {Evans}, J. and {Agrawal}, A. and {Leonard}, A.~J. and
    {Robitaille}, T.~P. and {Mampaey}, B. and {Iv{\'a}n Campos-Rozo}, J. and
    {Kirk}, M.~S.},
    title = "{SunPy{\mdash}Python for solar physics}",
  journal = {Computational Science and Discovery},
archivePrefix = "arXiv",
   eprint = {1505.02563},
 primaryClass = "astro-ph.IM",
     year = 2015,
    month = jan,
   volume = 8,
   number = 1,
      eid = {014009},
    pages = {014009},
      doi = {10.1088/1749-4699/8/1/014009},
   adsurl = {http://adsabs.harvard.edu/abs/2015CS%26D....8a4009S},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}"""


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

    from sunpy.util.logger import _init_log
    log = _init_log(config=config)

    __all__ = ['config', 'self_test', 'system_info']
