from functools import partial
import urllib2

import pytest

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
