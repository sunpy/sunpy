#!/usr/bin/env python
import os
import sys
import itertools

from setuptools.config import read_configuration

# Append cwd for pip 19
sys.path.append(os.path.abspath("."))
import ah_bootstrap  # noqa

from astropy_helpers.setup_helpers import setup, register_commands # noqa

################################################################################
# Override the default Astropy Test Command
################################################################################
cmdclass = register_commands()
try:
    from sunpy.tests.setup_command import SunPyTest
    # Overwrite the Astropy Testing framework
    cmdclass['test'] = type('SunPyTest', (SunPyTest,),
                            {'package_name': 'sunpy'})
except Exception:
    # Catch everything, if it doesn't work, we still want SunPy to install.
    pass

extras = read_configuration("setup.cfg")['options']['extras_require']
extras['all'] = list(itertools.chain.from_iterable(extras.values()))

setup(extras_require=extras)
