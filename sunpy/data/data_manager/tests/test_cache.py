from unittest.mock import patch

import pytest

from sunpy.util import hash_file
from sunpy.util.exceptions import SunpyUserWarning
from .mocks import MOCK_HASH


def test_cache_basic(cache):
    cache.download('https://example.com/abc.text')
    assert cache._downloader.times_called == 1


def test_cache_caching(cache):
    cache.download('https://example.com/abc.text')
    cache.download('https://example.com/abc.text')
    assert cache._downloader.times_called == 1


def test_cache_redownload(cache):
    cache.download('https://example.com/abc.text')
    cache.download('https://example.com/abc.text', redownload=True)
    assert cache._downloader.times_called == 2


def test_get_by_url(cache):
    cache.download('https://example.com/file_name')
    details = cache._get_by_url('https://example.com/file_name')
    assert details['file_path'].endswith('file_name')


def test_get_by_url_fail(cache):
    cache.download('https://example.com/file_name')
    details = cache._get_by_url('https://example.com/file')
    assert details is None


def test_get_by_hash(cache):
    cache.download('https://example.com/file_name')
    details = cache.get_by_hash(MOCK_HASH)
    assert details['file_path'].endswith('file_name')


def test_get_by_hash_fail(cache):
    cache.download('https://example.com/file_name')
    details = cache.get_by_hash('wrong_hash')
    assert details is None


def test_check_old_file_is_not_removed(cache, mocker):
    # Checks that https://github.com/sunpy/sunpy/issues/7249 is fixed

    # First check that the file is downloaded
    cache.download('https://example.com/file_name')
    first_details = cache.get_by_hash(MOCK_HASH)
    assert first_details['file_path'].endswith('file_name')

    # Force a redownload and check that the old file is not removed
    # We force an error inside `_download_and_hash` to test the full code path
    with patch('sunpy.data.data_manager.cache.hash_file') as download:
        download.side_effect = IOError
        with pytest.warns(SunpyUserWarning, match="Due to the above error, you"):
            path = cache.download('https://example.com/file_name', redownload=True)
        assert download.call_count == 1
    second_details = cache.get_by_hash(MOCK_HASH)
    assert first_details['file_path'] == second_details['file_path'] == str(path)


def test_file_change(cache, mocker):
    cache.download('https://example.com/abc.text')
    file_path = cache.get_by_hash(MOCK_HASH)["file_path"]
    initial_hash = hash_file(file_path)

    with open(file_path, "a") as file:
        file.write("test writing")

    modified_hash = hash_file(file_path)
    assert initial_hash != modified_hash

    with patch('sunpy.data.data_manager.cache.Cache._download_and_hash') as mock_download_and_hash:
        mock_download_and_hash.return_value = (file_path, initial_hash, 'https://example.com/abc.text')
        # Overwrite == True, so file should be overwritten with a new file
        file_path, file_hash, url = mock_download_and_hash('https://example.com/abc.text', redownload=True, overwrite=True)
        assert file_hash == initial_hash

        mock_download_and_hash.return_value = (file_path, modified_hash, "https://example.com/abc.text")
        #Overwrite == False, so the existing file should not be overwritten
        file_path, file_hash, url = mock_download_and_hash('https://example.com/abc.text', redownload=True)
        assert file_hash == modified_hash
