#!/usr/bin/env python
from setuptools import setup  # isort:skip
import os
import sys
from itertools import chain

from extension_helpers import get_extensions
from setuptools.config import read_configuration

################################################################################
# Raise helpful messages for test and build_docs commands
################################################################################
test_help = """\
Running tests is no longer done using 'python setup.py test'.

Instead you will need to run:
    tox -e offline
if you don't already have tox installed, you can install it with:
    pip install tox
if you only want to run part of the test suite, you can also use pytest directly with:
    pip install -e .[dev]
    pytest
for more information, see:
  https://docs.sunpy.org/en/latest/dev_guide/tests.html
"""

if 'test' in sys.argv:
    print(test_help)
    sys.exit(1)

docs_help = """\
Building the documentation is no longer done using 'python setup.py build_docs'.

Instead you will need to run:
    tox -e build_docs
if you don't already have tox installed, you can install it with:
    pip install tox
for more information, see:
   https://docs.sunpy.org/en/latest/dev_guide/documentation.html#usage
"""

if 'build_docs' in sys.argv or 'build_sphinx' in sys.argv:
    print(docs_help)
    sys.exit(1)

################################################################################
# Actual setup.py content
################################################################################

VERSION_TEMPLATE = """
# Note that we need to fall back to the hard-coded version if either
# setuptools_scm can't be imported or setuptools_scm can't determine the
# version, so we catch the generic 'Exception'.
try:
    from setuptools_scm import get_version
    __version__ = get_version(root='..', relative_to=__file__)
except Exception:
    __version__ = '{version}'
""".lstrip()

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

setup(extras_require=extras,
      use_scm_version={'write_to': os.path.join('sunpy', 'version.py'),
                       'write_to_template': VERSION_TEMPLATE},
      ext_modules=get_extensions())
