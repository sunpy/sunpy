from sunpy.data.manager.manager import DataManager
from sunpy.data.manager.storage import InMemStorage
from sunpy.data.manager.downloader import DownloaderBase
from sunpy.data.manager.cache import Cache

from pathlib import Path

import pytest


def write_to_test_file(contents):
    # TODO: Use tempfile. here and other places.
    with open('/tmp/test_file', 'w') as f:
        f.write(contents)


class MockDownloader(DownloaderBase):
    """
    MockDownloader.
    """

    def __init__(self):
        write_to_test_file("a")
        self.times_called = 0

    def download(self, url):
        self.times_called += 1
        return "/tmp/test_file"


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
    @manager.require('test_file', ['url1', 'url2'], '86f7e437faa5a7fce15d1ddcb9eaeaea377667b8')
    def foo(manager_tester=lambda x: 1):
        manager_tester(manager)

    return foo


def test_basic(manager, storage, downloader, data_function):
    data_function()

    assert downloader.times_called == 1
    assert len(storage._store) == 1
    assert storage._store[0]['file_path'] == '/tmp/test_file'


def test_cache(manager, storage, downloader, data_function):
    """
    Test calling function multiple times does not redownload.
    """
    data_function()
    data_function()

    assert downloader.times_called == 1
    assert len(storage._store) == 1
    assert storage._store[0]['file_path'] == '/tmp/test_file'


def test_skip_all(manager, storage, downloader, data_function):
    """
    Test skip_hash_check redownloads data.
    """
    data_function()
    with manager.skip_hash_check():
        data_function()

    assert downloader.times_called == 2
    assert len(storage._store) == 1
    assert storage._store[0]['file_path'] == '/tmp/test_file'


def test_replace_file(manager, storage, downloader, data_function):
    """
    Test the replace_file functionality.
    """

    def default_tester(manager):
        """
        Function to test whether the file is /tmp/test_file.
        """
        assert manager.get('test_file') == Path('/tmp/test_file')

    def replace_file_tester(manager):
        """
        Function to test whether the file is /tmp/lil.
        """
        assert manager.get('test_file') == Path('/tmp/lil')

    # Outside the context manager file is default
    data_function(default_tester)

    with manager.replace_file('test_file', 'file:///tmp/lil'):
        # Inside the file is replaced
        data_function(replace_file_tester)

    # Even after context manager call outside the file is default
    data_function(default_tester)


def test_wrong_hash_error(manager, storage):
    storage._store.append({
        'file_path': '/tmp/test_file',
        'file_hash': 'aa',
        'url': 'url1'
    })
    @manager.require('test_file', ['url1', 'url2'], 'asdf')
    def foo():
        pass
    with pytest.raises(KeyError):
        foo()


def test_file_changed(data_function):
    # Download the file first
    data_function()

    # The file was then locally changed
    write_to_test_file("asd")

    # Now it should error
    with pytest.raises(KeyError):
        data_function()
