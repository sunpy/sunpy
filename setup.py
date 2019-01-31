#!/usr/bin/env python
import os
import sys
import itertools

from setuptools.config import read_configuration

sys.path.append(os.path.abspath("."))
import ah_bootstrap  # noqa

from astropy_helpers.setup_helpers import setup # noqa

# Create an 'all' extra which is all the others combined
extras = read_configuration("setup.cfg")['options']['extras_require']
extras['all'] = list(itertools.chain.from_iterable(extras.values()))

setup(extras_require=extras)
