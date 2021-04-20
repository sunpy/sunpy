import json
import socket
from os.path import isfile
from urllib.error import HTTPError
from urllib.request import urlopen

import pytest
import requests


def pytest_generate_tests(metafunc):
    funcarglist = [None]
    if isfile(metafunc.config.getini("intercept_dump_file")):
        with open(metafunc.config.getini("intercept_dump_file")) as fd:
            funcarglist = json.load(fd)
            funcarglist = funcarglist.get(metafunc.function.__name__.replace("test_", ''), [None]) or [None]
    metafunc.parametrize("url", funcarglist, indirect=True)


@pytest.fixture
def url(request):
    if not request.param:
        pytest.skip("got empty parameter set")
    return request.param


@pytest.fixture(scope="function")
def skip_conditions(request):
    if not hasattr(request.config.option, "intercept_remote"):
        pytest.skip("pytest-intercept-remote plugin not loaded")

    if request.config.option.intercept_remote:
        pytest.skip("rerun without --intercept-remote option")

    if not isfile(request.config.getini("intercept_dump_file")):
        pytest.skip("intercept_dump_file not found")


@pytest.mark.remote_data
@pytest.mark.remote_status
@pytest.mark.usefixtures("skip_conditions")
def test_urls_urllib(url):
    try:
        res = urlopen(url)
        assert res.status == 200
    except HTTPError as e:
        pytest.xfail(f"URL unreachable, status:{e.code}")


@pytest.mark.remote_data
@pytest.mark.remote_status
@pytest.mark.usefixtures("skip_conditions")
def test_urls_requests(url):
    res = requests.get(url)
    status = res.status_code
    if status != 200:
        pytest.xfail(f"URL unreachable, status:{status}")
        return
    assert res.status_code == 200


@pytest.mark.remote_data
@pytest.mark.remote_status
@pytest.mark.usefixtures("skip_conditions")
def test_urls_socket(url):
    sock = socket.socket(socket.AF_INET)
    if len(url) == 4:
        sock = socket.socket(socket.AF_INET6)
    try:
        assert sock.connect(tuple(url)) is None
    except ConnectionRefusedError:
        pytest.xfail("URL unreachable")
    finally:
        sock.close()
