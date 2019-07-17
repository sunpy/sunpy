#!/usr/bin/env python
import os
import sys
from itertools import chain

from setuptools import setup
from setuptools.config import read_configuration

# Append cwd for pip 19
sys.path.append(os.path.abspath("."))
import ah_bootstrap  # noqa

from astropy_helpers.setup_helpers import register_commands, get_package_info # noqa

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

################################################################################
# Programmatically generate some extras combos.
################################################################################
extras = read_configuration("setup.cfg")['options']['extras_require']

# Dev is everything
extras['dev'] = list(chain(*extras.values()))

# All is everything but tests and docs
exclude_keys = ("tests", "docs", "dev")
ex_extras = dict(filter(lambda i: i[0] not in exclude_keys, extras.items()))
# Concatenate all the values together for 'all'
extras['all'] = list(chain.from_iterable(ex_extras.values()))

package_info = get_package_info()
setup(extras_require=extras, use_scm_version=True, cmdclass=cmdclass, **package_info)
