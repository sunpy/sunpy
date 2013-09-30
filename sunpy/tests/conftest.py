from functools import partial
import urllib2

import pytest


def site_reachable(url):
    try:
        urllib2.urlopen(url, timeout=1)
    except urllib2.URLError:
        return False
    except:
        raise
    else:
        return True


# Using a numerical IP-address avoids a DNS lookup, which may block the
# urllib2.urlopen call for more than a second
GOOGLE_URL = 'http://www.google.com'

is_online = partial(site_reachable, GOOGLE_URL)


def pytest_runtest_setup(item):
    """pytest hook to skip all tests that have the mark 'online' if the
    client is online (simply detected by checking whether http://www.google.com
    can be requested).

    """
    if isinstance(item, item.Function):
        if 'online' in item.keywords and not is_online():
            msg = 'skipping test {} (reason: client seems to be offline)'
            pytest.skip(msg.format(item.name))
