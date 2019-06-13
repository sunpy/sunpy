import pytest


def test_find_by_key_success(sqlstorage):
    test_details = {
        'file_hash': 'hash1',
        'file_path': '/tmp/test_file1',
        'url': 'http://example.com/test_file_1'
    }
    details = sqlstorage.find_by_key('file_hash', 'hash1')
    assert details == test_details

    details = sqlstorage.find_by_key('url', 'http://example.com/test_file_1')
    assert details == test_details


def test_find_by_key_fail(sqlstorage):
    details = sqlstorage.find_by_key('url', 'no_exist')
    assert details is None


def test_invalid_key(sqlstorage):
    with pytest.raises(KeyError):
        sqlstorage.find_by_key('key_not', 'hash_1')
