"""
SunPy
=====

An open-source Python library for Solar Physics data analysis.

* Homepage: https://sunpy.org
* Documentation: https://docs.sunpy.org/en/stable/
"""
import os
import sys
import logging

from sunpy.tests.runner import SunPyTestRunner
from sunpy.util import system_info
from sunpy.util.config import load_config, print_config
from sunpy.util.logger import _init_log
from .version import __version__

# Enforce Python version check during package import.
__minimum_python_version__ = "3.6"


class UnsupportedPythonError(Exception):
    pass


if sys.version_info < tuple(int(val) for val in __minimum_python_version__.split('.')):
    # This has to be .format to keep backwards compatibly.
    raise UnsupportedPythonError(
        "Sunpy does not support Python < {}".format(__minimum_python_version__))


def _get_bibtex():
    # Set the bibtex entry to the article referenced in CITATION.rst
    citation_file = os.path.join(os.path.dirname(__file__), 'CITATION.rst')

    with open(citation_file, 'r') as citation:
        refs = citation.read().split('@ARTICLE')[1:]
        if len(refs) == 0:
            return ''
        bibtexreference = r"@ARTICLE{}".format(refs[0])
    return bibtexreference


__citation__ = __bibtex__ = _get_bibtex()


self_test = SunPyTestRunner.make_test_runner_in(os.path.dirname(__file__))

# Load user configuration
config = load_config()


# Use the root logger as a dummy log before initializing Astropy's logger
log = logging.getLogger()

log = _init_log(config=config)

__all__ = ['config', 'self_test', 'system_info']
