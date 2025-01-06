import os
import re
from pathlib import Path

import pytest

from sunpy.data.data_manager.tests.mocks import MOCK_HASH, write_to_test_file
from sunpy.util.exceptions import SunpyUserWarning


def test_basic(storage, downloader, data_function):
    data_function()

    assert downloader.times_called == 1
    assert len(storage._store) == 1
    assert Path(storage._store[0]['file_path']).name == ('sunpy.test_file')


def test_download_cache(manager, storage, downloader, data_function):
    """
    Test calling function multiple times does not redownload.
    """
    data_function()
    data_function()

    assert downloader.times_called == 1
    assert len(storage._store) == 1
    assert Path(storage._store[0]['file_path']).name == ('sunpy.test_file')


def test_file_tampered(manager, storage, downloader, data_function):
    """
    Test calling function multiple times does not redownload.
    """
    data_function()
    write_to_test_file(manager._tempdir + '/sunpy.test_file', 'b')
    with pytest.warns(SunpyUserWarning):
        data_function()

    assert downloader.times_called == 2
    assert len(storage._store) == 1
    assert Path(storage._store[0]['file_path']).name == ('sunpy.test_file')


def test_wrong_hash_provided(manager):
    @manager.require('test_file', ['url1'], 'wrong_hash')
    def test_foo():
        pass

    with pytest.raises(RuntimeError):
        test_foo()


def test_defer_download(manager, storage, downloader, data_function, tmpdir):
    """
    Test that files are not downloaded immediately if defer_download is True,
    but are downloaded when get is called.
    """
    folder = tmpdir.strpath
    @manager.require('test_file', [f'file://{folder}/another_file'], MOCK_HASH,defer_download=True)
    def deferred_function():
        pass

    deferred_function()
    assert downloader.times_called == 0
    assert len(storage._store) == 0

def test_defer_download_get(manager, storage, downloader, data_function, tmpdir):
    folder = tmpdir.strpath
    @manager.require('test_file', [f'file://{folder}/another_file'], MOCK_HASH, defer_download=True)
    def deferred_function():
        manager.get('test_file')

    deferred_function()
    assert downloader.times_called == 1
    assert len(storage._store) == 1


def test_skip_all(manager, storage, downloader, data_function):
    """
    Test skip_hash_check redownloads data.
    """
    data_function()
    with manager.skip_hash_check():
        data_function()

    assert downloader.times_called == 2
    assert len(storage._store) == 1
    assert Path(storage._store[0]['file_path']).name == ('sunpy.test_file')


def test_override_file(manager, storage, downloader, data_function, tmpdir):
    """
    Test the override_file functionality.
    """

    def default_tester(manager):
        """
        Function to test whether the file name is test_file.
        """
        assert manager.get('test_file').name == ('sunpy.test_file')

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

    with pytest.raises(ValueError, match="Hash provided to override_file does not match hash of the file."):
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
    with pytest.raises(ValueError, match=re.escape("['url1', 'url2'] has already been downloaded, but no file matching the hash asdf can be found.")):
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


def test_delete_db(sqlmanager, sqlstorage):
    # Download the file
    @sqlmanager.require('test_file', ['http://example.com/test_file'], MOCK_HASH)
    def test_function():
        pass

    test_function()

    # The DB file was then deleted
    os.remove(str(sqlstorage._db_path))

    # SQLite should not throw an error
    test_function()


def test_same_file_id_different_module(downloader, storage,
                                       data_function, data_function_from_fake_module):
    # Uses name 'test_file' to refer to the file
    data_function()

    # Change hash of the above file to allow MockDownloader to download another file
    # Otherwise it will skip the download because a file with the same hash already exists
    storage._store[0]['file_hash'] = 'abc'

    # This function from a different module uses same name 'test_file' to refer to a different file
    data_function_from_fake_module()

    assert len(storage._store) == 2
    assert downloader.times_called == 2

    # Check if the files are namespaced correctly
    assert Path(storage._store[0]['file_path']).name == 'sunpy.test_file'
    assert Path(storage._store[1]['file_path']).name == 'fake_module.test_file'


def test_namespacing_with_manager_override_file(module_patched_manager, downloader,
                                                storage, data_function_from_fake_module):
    # Download a file using manager.require()
    data_function_from_fake_module()

    assert len(storage._store) == 1
    assert downloader.times_called == 1
    assert Path(storage._store[0]['file_path']).name == 'fake_module.test_file'

    # Override the file name with a different URI
    with module_patched_manager.override_file(
            'test_file', 'http://www.different_uri.com/new_file', MOCK_HASH):
        data_function_from_fake_module()

        assert downloader.times_called == 2

        # New file entry is stored in manager._file_cache only
        # It's not stored in InMemStorage or SqlStorage
        assert len(storage._store) == 1

        assert Path(
            module_patched_manager._file_cache['test_file']['fake_module.']
        ).name == 'fake_module.new_file'

        # Storage still contains original test_file
        assert Path(storage._store[0]['file_path']).name == 'fake_module.test_file'

    # Request the original file again
    data_function_from_fake_module()

    # File doesn't get redownloaded, instead it is retrieved using the file hash
    assert downloader.times_called == 2

    # new_file entry in manager._file_cache is replaced with the original test_file
    assert Path(
        module_patched_manager._file_cache['test_file']['fake_module.']
    ).name == 'fake_module.test_file'

    # Storage still contains original test_file
    assert Path(storage._store[0]['file_path']).name == 'fake_module.test_file'


def test_file_deleted_redownload(storage, downloader, data_function):
    # Checks that if a file is deleted, a warning is raised and the file is re-downloaded.

    data_function()
    assert len(storage._store) == 1
    test_file_path = Path(storage._store[0]['file_path'])
    assert test_file_path.exists()

    os.remove(test_file_path)
    assert not test_file_path.exists()

    with pytest.warns(SunpyUserWarning, match="Requested file appears to missing and will be redownloaded."):
        data_function()

    assert downloader.times_called == 2
    assert len(storage._store) == 1
    assert test_file_path.exists()
