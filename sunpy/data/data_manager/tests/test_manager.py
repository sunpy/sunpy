from pathlib import Path

import pytest

from sunpy.data.data_manager.tests.mocks import MOCK_HASH, write_to_test_file
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


def test_override_file(manager, storage, downloader, data_function, tmpdir):
    """
    Test the override_file functionality.
    """

    def default_tester(manager):
        """
        Function to test whether the file name is test_file.
        """
        assert manager.get('test_file').name == ('test_file')

    def override_file_tester(manager):
        """
        Function to test whether the file is /tmp/another_file.
        """
        assert manager.get('test_file') == Path(f'{folder}/another_file')

    # Outside the context manager file is default
    folder = tmpdir.strpath
    data_function(default_tester)
    write_to_test_file(str(Path(folder+'/another_file')), 'a')

    with manager.override_file('test_file', f'file://{folder}/another_file'):
        # Inside the file is replaced
        data_function(override_file_tester)

    # TODO: this combined with the check above fails on windows
    # with manager.override_file('test_file', f'{folder}/another_file'):
    #     # Inside the file is replaced
    #     data_function(override_file_tester)

    # check the function works with hash provided
    with manager.override_file('test_file', f'file://{folder}/another_file', MOCK_HASH):
        data_function(override_file_tester)

    with pytest.raises(ValueError):
        # check if functions errors with the wrong hash
        with manager.override_file('test_file', f'file://{folder}/another_file', 'wrong_hash'):
            # Inside the file is replaced
            data_function(override_file_tester)

    # Even after context manager call outside the file is default
    data_function(default_tester)


def test_override_file_remote(manager, downloader, data_function):
    replace_url = 'http://example.com/another_file'
    data_function()
    assert downloader.times_called == 1
    with manager.override_file('test_file', replace_url):
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
    with pytest.raises(ValueError):
        foo()


def test_file_changed(data_function, storage):
    # Download the file first
    data_function()

    file = storage._store[0]['file_path']

    # The file was then locally changed
    write_to_test_file(file, "asd")

    # Now it should error
    with pytest.warns(SunpyUserWarning):
        data_function()
