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
from sunpy.tests.helpers import new_hash_library, figure_test_pngfiles
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
    if len(figure_test_pngfiles) > 0:
        tempdir = tempfile.mkdtemp(suffix="_figures")

        # Rename each PNG with the name of the corresponding test
        for test_name in figure_test_pngfiles:
            os.rename(figure_test_pngfiles[test_name], os.path.join(tempdir, test_name + '.png'))

        # Write the new hash library in JSON
        hashfile = os.path.join(tempdir, HASH_LIBRARY_NAME)
        with open(hashfile, 'w') as outfile:
            json.dump(new_hash_library, outfile, sort_keys=True, indent=4, separators=(',', ': '))

        print('All test files for figure hashes can be found in {0}'.format(tempdir))
        print("The corresponding hash library is {0}".format(hashfile))
