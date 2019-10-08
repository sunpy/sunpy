from pathlib import Path

import pytest

from sunpy.data.data_manager.tests.mocks import write_to_test_file, MOCK_HASH
from sunpy.util.exceptions import SunpyUserWarning


def test_basic(storage, downloader, data_function):
    data_function()

    assert downloader.times_called == 1
    assert len(storage._store) == 1
    assert Path(storage._store[0]['file_path']).name == ('test_file')


def test_cache(manager, storage, downloader, data_function):
    """
    Test calling function multiple times does not redownload.
    """
    data_function()
    data_function()

    assert downloader.times_called == 1
    assert len(storage._store) == 1
    assert Path(storage._store[0]['file_path']).name == ('test_file')


def test_file_tampered(manager, storage, downloader, data_function):
    """
    Test calling function multiple times does not redownload.
    """
    data_function()
    write_to_test_file(manager._tempdir + '/test_file', 'b')
    with pytest.warns(SunpyUserWarning):
        data_function()

    assert downloader.times_called == 2
    assert len(storage._store) == 1
    assert Path(storage._store[0]['file_path']).name == ('test_file')


def test_wrong_hash_provided(manager):
    @manager.require('test_file', ['url1'], 'wrong_hash')
    def test_foo():
        pass

    with pytest.raises(RuntimeError):
        test_foo()

def test_skip_all(manager, storage, downloader, data_function):
    """
    Test skip_hash_check redownloads data.
    """
    data_function()
    with manager.skip_hash_check():
        data_function()

    assert downloader.times_called == 2
    assert len(storage._store) == 1
    assert Path(storage._store[0]['file_path']).name == ('test_file')

def test_replace_file(manager, storage, downloader, data_function):
    """
    Test the replace_file functionality.
    """

    def default_tester(manager):
        """
        Function to test whether the file name is test_file.
        """
        assert manager.get('test_file').name == ('test_file')

    def replace_file_tester(manager):
        """
        Function to test whether the file is /tmp/another_file.
        """
        assert manager.get('test_file') == Path('/tmp/another_fiel')

    # Outside the context manager file is default
    data_function(default_tester)
    write_to_test_file('/tmp/another_fiel', 'a')

    with manager.replace_file('test_file', 'file:///tmp/another_fiel'):
        # Inside the file is replaced
        data_function(replace_file_tester)

    # check the function works with hash provided
    with manager.replace_file('test_file', 'file:///tmp/another_fiel', MOCK_HASH):
        data_function(replace_file_tester)

    with pytest.raises(KeyError):
        # check if functions errors with the wrong hash
        with manager.replace_file('test_file', 'file:///tmp/another_fiel', 'wrong_hash'):
            # Inside the file is replaced
            data_function(replace_file_tester)

    # Even after context manager call outside the file is default
    data_function(default_tester)


def test_replace_file_remote(manager, downloader, data_function):
    replace_url = 'http://example.com/another_fiel'
    data_function()
    assert downloader.times_called == 1
    with manager.replace_file('test_file', replace_url):
        data_function()

    assert downloader.times_called == 2
    assert downloader.last_called_url == replace_url


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


def test_file_changed(data_function, storage):
    # Download the file first
    data_function()

    file = storage._store[0]['file_path']

    # The file was then locally changed
    write_to_test_file(file, "asd")

    # Now it should error
    with pytest.raises(KeyError):
        data_function()
