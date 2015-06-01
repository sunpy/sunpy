from functools import partial
import urllib2
import os
import tempfile
import json

import pytest

# Force MPL to use non-gui backends for testing.
try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

from sunpy.tests import hash

hash_library_original_len = len(hash.hash_library)

GOOGLE_URL = 'http://www.google.com'


def site_reachable(url):
    try:
        urllib2.urlopen(url, timeout=1)
    except urllib2.URLError:
        return False
    else:
        return True


is_online = partial(site_reachable, GOOGLE_URL)


def pytest_runtest_setup(item):
    """pytest hook to skip all tests that have the mark 'online' if the
    client is online (simply detected by checking whether http://www.google.com
    can be requested).

    """
    if isinstance(item, item.Function):
        if 'online' in item.keywords and not is_online():
            msg = 'skipping test {0} (reason: client seems to be offline)'
            pytest.skip(msg.format(item.name))

def pytest_unconfigure(config):
    #Check if additions have been made to the hash library
    if len(hash.hash_library) > hash_library_original_len:
        #Write the new hash library in JSON
        tempdir = tempfile.mkdtemp()
        hashfile = os.path.join(tempdir, hash.HASH_LIBRARY_NAME)
        with open(hashfile, 'wb') as outfile:
            json.dump(hash.hash_library, outfile, sort_keys=True, indent=4, separators=(',', ': '))
        print("The hash library has expanded and should be copied to sunpy/tests/")
        print("  " + hashfile)
