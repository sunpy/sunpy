from __future__ import absolute_import, print_function
from functools import partial

import os
import tempfile
import json

# Force MPL to use non-gui backends for testing.
try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

from sunpy.tests.hash import HASH_LIBRARY_NAME
from sunpy.tests.helpers import new_hash_library, test_fig_dir
from sunpy.extern import six

import pytest


# Don't actually import pytest_remotedata because that can do things to the
# entrypoints code in pytest.
if six.PY2:
    import imp
    try:
        imp.find_module('pytest_remotedata')
        HAVE_REMOTEDATA = True
    except ImportError:
        HAVE_REMOTEDATA = False
else:
    import importlib
    remotedata_spec = importlib.util.find_spec("pytest_remotedata")
    HAVE_REMOTEDATA = remotedata_spec is not None


def pytest_runtest_setup(item):
    """
    pytest hook to skip all tests that have the mark 'online' if the
    client is online (simply detected by checking whether http://www.google.com
    can be requested).
    """
    if isinstance(item, item.Function):
        if 'remote_data' in item.keywords and not HAVE_REMOTEDATA:
            pytest.skip("skipping remotedata tests as pytest-remotedata is not installed")


def pytest_unconfigure(config):
    if len(new_hash_library) > 0:
        # Write the new hash library in JSON
        hashfile = os.path.join(test_fig_dir, HASH_LIBRARY_NAME)
        with open(hashfile, 'w') as outfile:
            json.dump(new_hash_library, outfile, sort_keys=True, indent=4, separators=(',', ': '))

        print('All images from image tests can be found in {0}'.format(test_fig_dir))
        print("The corresponding hash library is {0}".format(hashfile))
