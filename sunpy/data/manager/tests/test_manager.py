from sunpy.data.manager.manager import DataManager
from sunpy.data.manager.storage import InMemStorage
from sunpy.data.manager.downloader import DownloaderBase
from sunpy.data.manager.cache import Cache

from pathlib import Path

import pytest


class MockDownloader(DownloaderBase):
    """MockDownloader"""

    def __init__(self):
        self.times_called = 0

    def download(self, url):
        self.times_called += 1
        return "/tmp/lol"


@pytest.fixture
def downloader():
    downloader = MockDownloader()
    return downloader


@pytest.fixture
def storage():
    storage = InMemStorage()
    return storage


@pytest.fixture
def manager(downloader, storage):
    manager = DataManager(Cache(downloader, storage))
    return manager


@pytest.fixture
def data_function(manager):
    @manager.require('test_file', ['url1', 'url2'], 'asdf')
    def foo(manager_tester=lambda x: 1):
        manager_tester(manager)

    return foo


def test_basic(manager, storage, downloader, data_function):
    data_function()

    assert downloader.times_called == 1
    assert len(storage._store) == 1
    assert storage._store[0]['file_path'] == '/tmp/lol'


def test_cache(manager, storage, downloader, data_function):
    """Test calling function multiple times does not redownload"""
    data_function()
    data_function()

    assert downloader.times_called == 1
    assert len(storage._store) == 1
    assert storage._store[0]['file_path'] == '/tmp/lol'


def test_skip_all(manager, storage, downloader, data_function):
    """Test skip_hash_check redownloads data"""
    data_function()
    with manager.skip_hash_check():
        data_function()

    assert downloader.times_called == 2
    assert len(storage._store) == 1
    assert storage._store[0]['file_path'] == '/tmp/lol'


def test_replace_file(manager, storage, downloader, data_function):
    """Test the replace_file functionality"""

    def default_tester(manager):
        """Function to test whether the file is /tmp/lol"""
        assert manager.get('test_file') == Path('/tmp/lol')

    def replace_file_tester(manager):
        """Function to test whether the file is /tmp/lil"""
        assert manager.get('test_file') == Path('/tmp/lil')

    # Outside the context manager file is default
    data_function(default_tester)

    with manager.replace_file('test_file', 'file:///tmp/lil'):
        # Inside the file is replaced
        data_function(replace_file_tester)

    # Even after context manager call outside the file is default
    data_function(default_tester)
