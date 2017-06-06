from __future__ import absolute_import, print_function
from functools import partial

import os
import socket
import tempfile
import json

from sunpy.extern.six.moves.urllib.request import urlopen
from sunpy.extern.six.moves.urllib.error import URLError

import pytest


# Force MPL to use non-gui backends for testing.
try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

from astropy.tests import disable_internet

from sunpy.tests.hash import HASH_LIBRARY_NAME
from sunpy.tests.helpers import new_hash_library, figure_test_pngfiles

GOOGLE_URL = 'http://www.google.com'


def site_reachable(url):
    try:
        urlopen(url, timeout=1)
    except (URLError, socket.timeout):
        return False
    else:
        return True


is_online = partial(site_reachable, GOOGLE_URL)


def pytest_runtest_setup(item):
    """
    pytest hook to skip all tests that have the mark 'online' if the
    client is online (simply detected by checking whether http://www.google.com
    can be requested).
    """
    if isinstance(item, item.Function):
        if 'online' in item.keywords and not is_online():
            msg = 'skipping test {0} (reason: client seems to be offline)'
            pytest.skip(msg.format(item.name))

        if 'online' not in item.keywords:
            disable_internet.turn_off_internet()


def pytest_runtest_teardown(item, nextitem):
    disable_internet.turn_on_internet()


def pytest_unconfigure():
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
