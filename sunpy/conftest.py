from __future__ import absolute_import, print_function
from functools import partial

import os
import pathlib
import tempfile
import json

# Force MPL to use non-gui backends for testing.
try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

import sunpy.tests.helpers
from sunpy.tests.hash import HASH_LIBRARY_NAME
from sunpy.tests.helpers import new_hash_library, generate_figure_webpage
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


def pytest_addoption(parser):
    parser.addoption("--figure_dir", action="store", default="./figure_test_images")


@pytest.fixture(scope='session', autouse=True)
def figure_base_dir(request):
    sunpy.tests.helpers.figure_base_dir = pathlib.Path(
        request.config.getoption("--figure_dir"))


def pytest_runtest_setup(item):
    """
    pytest hook to skip all tests that have the mark 'remotedata' if the
    pytest_remotedata plugin is not installed.
    """
    if isinstance(item, item.Function):
        if 'remote_data' in item.keywords and not HAVE_REMOTEDATA:
            pytest.skip("skipping remotedata tests as pytest-remotedata is not installed")


def pytest_unconfigure(config):

    # If at least one figure test has been run, print result image directory
    if len(new_hash_library) > 0:
        # Write the new hash library in JSON
        figure_base_dir = pathlib.Path(config.getoption("--figure_dir"))
        hashfile = figure_base_dir / HASH_LIBRARY_NAME
        with open(hashfile, 'w') as outfile:
            json.dump(new_hash_library, outfile, sort_keys=True, indent=4, separators=(',', ': '))

        """
        Turn on internet when generating the figure comparison webpage.
        """
        if HAVE_REMOTEDATA:
            from pytest_remotedata.disable_internet import turn_on_internet, turn_off_internet
        else:
            def turn_on_internet(): pass
            def turn_off_internet(): pass

        turn_on_internet()
        generate_figure_webpage(new_hash_library)
        turn_off_internet()

        print('All images from image tests can be found in {0}'.format(figure_base_dir))
        print("The corresponding hash library is {0}".format(hashfile))
