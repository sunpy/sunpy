import os
import json
import pathlib
import warnings
import importlib
import tempfile

import pytest

import astropy
import astropy.config.paths

import sunpy.tests.helpers
from sunpy.tests.hash import HASH_LIBRARY_NAME
from sunpy.tests.helpers import new_hash_library, generate_figure_webpage
from sunpy.util.exceptions import SunpyDeprecationWarning

# Force MPL to use non-gui backends for testing.
try:
    import matplotlib
except ImportError:
    pass
else:
    matplotlib.use('Agg')

# Don't actually import pytest_remotedata because that can do things to the
# entrypoints code in pytest.
remotedata_spec = importlib.util.find_spec("pytest_remotedata")
HAVE_REMOTEDATA = remotedata_spec is not None


def pytest_addoption(parser):
    parser.addoption("--figure_dir", action="store", default="./figure_test_images")


@pytest.fixture(scope='session', autouse=True)
def figure_base_dir(request):
    sunpy.tests.helpers.figure_base_dir = pathlib.Path(
        request.config.getoption("--figure_dir"))


@pytest.fixture(scope='session', autouse=True)
def tmp_config_dir(request):
    """
    Globally set the default config for all tests.
    """
    tmpdir = tempfile.TemporaryDirectory()

    os.environ["SUNPY_CONFIGDIR"] = str(tmpdir.name)
    astropy.config.paths.set_temp_config._temp_path = str(tmpdir.name)
    astropy.config.paths.set_temp_cache._temp_path = str(tmpdir.name)

    yield

    del os.environ["SUNPY_CONFIGDIR"]
    astropy.config.paths.set_temp_config._temp_path = None
    astropy.config.paths.set_temp_cache._temp_path = None


@pytest.fixture()
def sunpy_cache(mocker):
    """
    Provide a way to add local files to the cache. This can be useful when mocking
    remote requests.
    """
    from types import MethodType
    from sunpy.data.data_manager.cache import Cache
    from sunpy.data.data_manager.storage import InMemStorage
    from sunpy.data.data_manager.downloader import ParfiveDownloader
    cache_dir = tempfile.mkdtemp()
    cache = Cache(
        ParfiveDownloader(),
        InMemStorage(),
        cache_dir,
        None
    )

    def add(self, url, path):
        self._storage.store({
            'url': url,
            'file_path': path,
            'file_hash': 'none',  # hash doesn't matter
        })
    cache.add = MethodType(add, cache)

    def func(mocked):
        mocker.patch(mocked, cache)
        return cache
    yield func


@pytest.fixture()
def undo_config_dir_patch():
    """
    Provide a way for certain tests to not have the config dir.
    """
    oridir = os.environ["SUNPY_CONFIGDIR"]
    del os.environ["SUNPY_CONFIGDIR"]
    yield
    os.environ["SUNPY_CONFIGDIR"] = oridir


@pytest.fixture(scope='session', autouse=True)
def tmp_dl_dir(request):
    """
    Globally set the default download directory for the test run to a tmp dir.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        os.environ["SUNPY_DOWNLOADDIR"] = tmpdir
        yield tmpdir
        del os.environ["SUNPY_DOWNLOADDIR"]


@pytest.fixture()
def undo_download_dir_patch():
    """
    Provide a way for certain tests to not have tmp download dir.
    """
    oridir = os.environ["SUNPY_DOWNLOADDIR"]
    del os.environ["SUNPY_DOWNLOADDIR"]
    yield
    os.environ["SUNPY_DOWNLOADDIR"] = oridir


def pytest_runtest_setup(item):
    """
    pytest hook to skip all tests that have the mark 'remotedata' if the
    pytest_remotedata plugin is not installed.
    """
    if isinstance(item, pytest.Function):
        if 'remote_data' in item.keywords and not HAVE_REMOTEDATA:
            pytest.skip("skipping remotedata tests as pytest-remotedata is not installed")


def pytest_configure(config):
    """
    Used to detect if code is being run as a test or not
    From: https://docs.pytest.org/en/5.1.1/example/simple.html#detect-if-running-from-within-a-pytest-run
    """
    import sunpy
    sunpy._called_from_test = True


def pytest_unconfigure(config):
    # tear down for setup in pytest configure
    import sunpy
    del sunpy._called_from_test

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

        print('All images from image tests can be found in {}'.format(figure_base_dir.resolve()))
        print("The corresponding hash library is {}".format(hashfile.resolve()))


def pytest_sessionstart(session):
    warnings.simplefilter("error", SunpyDeprecationWarning)
