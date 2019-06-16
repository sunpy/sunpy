from pathlib import Path

import pytest

from sunpy.data.data_manager.tests.mocks import write_to_test_file


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
        Function to test whether the file is /tmp/test_file.
        """
        assert manager.get('test_file').name == ('test_file')

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


def test_file_changed(data_function, storage):
    # Download the file first
    data_function()

    file = storage._store[0]['file_path']

    # The file was then locally changed
    write_to_test_file(file, "asd")

    # Now it should error
    with pytest.raises(KeyError):
        data_function()
