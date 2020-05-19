import pytest


def test_find_by_key_success(sqlstorage):
    test_details = {
        'file_hash': 'hash1',
        'file_path': '/tmp/test_file1',
        'url': 'http://example.com/test_file_1',
        'time': '2019-06-17T19:16:55.159274',
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


def test_insert(sqlstorage):
    details = {
        'file_hash': 'hash9',
        'file_path': '/tmp/test_file9',
        'url': 'http://example.com/test_file_9',
        'time': '2019-06-17T19:16:55.159274',
    }
    sqlstorage.store(details)
    new_details = sqlstorage.find_by_key('file_hash', 'hash9')
    assert details == new_details


def test_delete(sqlstorage):
    sqlstorage.delete_by_key('file_hash', 'hash1')
    details = sqlstorage.find_by_key('file_hash', 'hash1')
    assert details is None
