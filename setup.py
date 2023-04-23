#!/usr/bin/env python
from setuptools import setup  # isort:skip
import os

from extension_helpers import get_extensions

setup(
    use_scm_version={'write_to': os.path.join('sunpy', '_version.py')},
    ext_modules=get_extensions(),
)
