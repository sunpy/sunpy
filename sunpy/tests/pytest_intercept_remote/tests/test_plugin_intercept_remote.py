import socket
from urllib.request import urlopen

import pytest
import requests


@pytest.fixture(scope="function")
def skip_condition(request):
    if not hasattr(request.config.option, "intercept_remote"):
        pytest.skip("intercept-remote plugin not loaded")
    if request.config.option.intercept_remote:
        pytest.skip("rerun without --intercept-remote option")


@pytest.mark.remote_data
@pytest.mark.usefixtures("skip_condition")
def test_requests_urls():
    u = requests.get("https://www.python.org")
    assert u.status_code == 200


@pytest.mark.remote_data
@pytest.mark.usefixtures("skip_condition")
def test_urllib_urls():
    u = urlopen("https://www.python.org/")
    assert u.status == 200


@pytest.mark.remote_data
@pytest.mark.usefixtures("skip_condition")
def test_socket():
    s = socket.socket()
    assert s.connect(("www.python.org", 80)) is None
    s.close()


@pytest.mark.remote_data
@pytest.mark.usefixtures("skip_condition")
def test_remote(testdir):
    testdir.makeconftest(
        """
        pytest_plugins=["sunpy.tests.pytest_intercept_remote.plugin"]
        """)
    testdir.makepyfile(
        """
        import pytest
        from urllib.request import urlopen
        import requests
        import socket

        @pytest.mark.remote_data
        def test_requests_urls():
            u = requests.get("https://www.python.org")
            assert u.status_code == 200

        @pytest.mark.remote_data
        def test_urllib_urls():
            u = urlopen("https://www.python.org/")
            assert u.status == 200

        @pytest.mark.remote_data
        def test_socket():
            s = socket.socket()
            assert s.connect(("www.python.org", 80)) is None
            s.close()

        @pytest.mark.remote_data
        def test_dump(intercepted_urls):
            assert intercepted_urls == {"urls_urllib": [], "urls_requests": [], "urls_socket": []}
        """
    )

    result = testdir.runpytest("-q", "-p", "no:warnings", "--remote-data=any",
                               "-o", "intercept_dump_file=test_urls.json")
    result.assert_outcomes(passed=4)


@pytest.mark.remote_data
@pytest.mark.usefixtures("skip_condition")
def test_intercept_remote(testdir):
    testdir.makeconftest(
        """
        pytest_plugins=["sunpy.tests.pytest_intercept_remote.plugin"]
        """)
    testdir.makepyfile(
        """
        import pytest
        from urllib.request import urlopen
        import requests
        import socket

        @pytest.mark.remote_data
        def test_requests_urls():
            u = requests.get("https://www.python.org")
            assert u.status_code == 200

        @pytest.mark.remote_data
        def test_urllib_urls():
            u = urlopen("https://www.python.org/")
            assert u.status == 200

        @pytest.mark.remote_data
        def test_socket():
            s = socket.socket()
            assert s.connect(("www.python.org", 80)) is None
            s.close()

        @pytest.mark.remote_data
        def test_dump(intercepted_urls):
            assert intercepted_urls == {"urls_urllib": ["https://www.python.org/"],
                                        "urls_requests": ["https://www.python.org/"],
                                        "urls_socket": [("www.python.org", 80)]}
        """
    )

    result = testdir.runpytest("-q", "-p", "no:warnings", "--remote-data=any", "--intercept-remote",
                               "-o", "intercept_dump_file=test_urls.json")
    result.assert_outcomes(xfailed=3, passed=1)
