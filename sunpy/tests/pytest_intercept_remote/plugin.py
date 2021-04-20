import json

import pytest

mpatch = pytest.MonkeyPatch()

_requests_urls = []
_sockets_urls = []
_urllib_urls = []


def pytest_addoption(parser):
    DEFAULT_DUMP_FILE = "sunpy/tests/remote_urls.json"

    parser.addoption("--intercept-remote", action="store_true", default=False,
                     help="Intercepts outgoing connections requests.")
    parser.addini("intercept_dump_file", "filepath at which intercepted requests are dumped",
                  type="string", default=DEFAULT_DUMP_FILE)


def pytest_configure(config):
    if not config.option.intercept_remote and config.option.verbose:
        print("Intercept outgoing requests: disabled")

    remote_data = config.getoption("remote_data")
    intercept_remote = config.getoption('--intercept-remote')

    if remote_data and intercept_remote:
        global mpatch
        intercept_patch(mpatch)


def pytest_unconfigure(config):
    """
    Dump requests and clean
    """
    if config.option.intercept_remote:
        global mpatch
        mpatch.undo()
        intercept_dump(config)


def urlopen_mock(self, http_class, req, **http_conn_args):
    """
    Mock function for urllib.request.urlopen.
    """
    global _urllib_urls
    _urllib_urls.append(req.get_full_url())
    pytest.xfail(f"The test was about to call {req.get_full_url()}")


def requests_mock(self, method, url, *args, **kwargs):
    """
    Mock function for urllib3 module.
    """
    global _requests_urls
    full_url = f"{self.scheme}://{self.host}{url}"
    _requests_urls.append(full_url)
    pytest.xfail(f"The test was about to {method} {full_url}")


def socket_connect_mock(self, addr):
    """
    Mock function for socket.socket.
    """
    global _sockets_urls
    self.close()
    host = addr[0]
    port = addr[1]
    _sockets_urls.append(addr)
    pytest.xfail(f"The test was about to connect to {host}:{port}")


def intercept_patch(mpatch):
    """
    Monkey Patches urllib, urllib3 and socket.
    """
    mpatch.setattr(
        "urllib.request.AbstractHTTPHandler.do_open", urlopen_mock)
    mpatch.setattr(
        "urllib3.connectionpool.HTTPConnectionPool.urlopen", requests_mock)
    mpatch.setattr(
        "socket.socket.connect", socket_connect_mock)


@pytest.fixture
def intercepted_urls():
    """
    Pytest fixture for getting intercepted urls in a test
    """
    _urls = {
        'urls_urllib': _urllib_urls,
        'urls_requests': _requests_urls,
        'urls_socket': _sockets_urls}
    return _urls


def intercept_dump(config):
    """
    Dumps intercepted requests to ini option ``intercept_dump_file``.
    """
    global _requests_urls, _urllib_urls, _sockets_urls

    _urls = {
        'urls_urllib': _urllib_urls,
        'urls_requests': _requests_urls,
        'urls_socket': _sockets_urls}
    with open(config.getini("intercept_dump_file"), 'w') as fd:
        json.dump(_urls, fd)
