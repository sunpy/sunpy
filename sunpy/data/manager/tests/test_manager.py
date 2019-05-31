from sunpy.data.manager.manager import DataManager
from sunpy.data.manager.storage import InMemStorage
from sunpy.data.manager.downloader import DownloaderBase

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
    manager = DataManager(downloader, storage)
    return manager


def test_basic(manager, storage, downloader):
    @manager.require('test_file', ['url1', 'url2'], 'hash')
    def foo():
        assert manager.get('test_file') == '/tmp/lol'
    foo()

    assert downloader.times_called == 1
    assert len(storage._store) == 1
    assert storage._store[0]['file_path'] == '/tmp/lol'


def test_cache(manager, storage, downloader):
    @manager.require('test_file', ['url1', 'url2'], 'hash')
    def foo():
        assert manager.get('test_file') == '/tmp/lol'
    foo()
    foo()

    assert downloader.times_called == 1
    assert len(storage._store) == 1
    assert storage._store[0]['file_path'] == '/tmp/lol'


def test_skip_all(manager, storage, downloader):
    @manager.require('test_file', ['url1', 'url2'], 'hash')
    def foo():
        assert manager.get('test_file') == '/tmp/lol'
    foo()
    with manager.skip_hash_check():
        foo()

    assert downloader.times_called == 2
    assert len(storage._store) == 1
    assert storage._store[0]['file_path'] == '/tmp/lol'
