#!/usr/bin/env python
from itertools import chain

from setuptools import setup
from setuptools.config import read_configuration

from extension_helpers.setup_helpers import get_extensions

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

setup(extras_require=extras, use_scm_version=True, ext_modules=get_extensions())
