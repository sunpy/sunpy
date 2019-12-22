import csv
import tempfile
from pathlib import Path
from unittest import mock

import pytest

from sunpy.data.data_manager.cache import Cache
from sunpy.data.data_manager.manager import DataManager
from sunpy.data.data_manager.storage import InMemStorage, SqliteStorage
from sunpy.data.data_manager.tests import mocks

DB_TESTDATA_FILE = Path(__file__).parent / 'db_testdata.csv'


@pytest.fixture
def downloader():
    downloader = mocks.MockDownloader()
    return downloader


@pytest.fixture
def storage():
    storage = InMemStorage()
    return storage


@pytest.fixture
def sqlstorage():
    temp = tempfile.mktemp(suffix='.db')
    storage = SqliteStorage(temp)
    with open(DB_TESTDATA_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            storage.store(row)
    return storage


@pytest.fixture
def cache(tmp_path, downloader, storage, mocker):
    m = mock.Mock()
    m.headers = {'Content-Disposition': 'test_file'}
    mocker.patch('sunpy.data.data_manager.cache.urlopen', return_value=m)
    cache = Cache(downloader, storage, tmp_path)
    yield cache


@pytest.fixture
def manager(tmp_path, downloader, storage, mocker):
    manager = DataManager(Cache(downloader, storage, tmp_path))
    manager._tempdir = str(tmp_path)
    m = mock.Mock()
    m.headers = {'Content-Disposition': 'test_file'}
    mocker.patch('sunpy.data.data_manager.cache.urlopen', return_value=m)
    yield manager


@pytest.fixture
def data_function(manager):
    @manager.require('test_file', ['url1/test_file', 'url2'], mocks.MOCK_HASH)
    def foo(manager_tester=lambda x: 1):
        manager_tester(manager)

    return foo
