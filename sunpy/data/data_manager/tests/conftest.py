import shutil
import tempfile
import csv
from unittest import mock

import pytest
from sunpy.data.data_manager.tests import mocks

from sunpy.data.data_manager.cache import Cache
from sunpy.data.data_manager.manager import DataManager
from sunpy.data.data_manager.storage import InMemStorage, SqliteStorage

DB_TESTDATA_FILE = 'sunpy/data/data_manager/tests/db_testdata.csv'


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
def cache(downloader, storage, mocker):
    tempdir = tempfile.mkdtemp()
    m = mock.Mock()
    m.headers = {'Content-Disposition': 'test_file'}
    mocker.patch('sunpy.data.data_manager.cache.urlopen', return_value=m)
    cache = Cache(downloader, storage, tempdir)
    yield cache
    shutil.rmtree(tempdir)


@pytest.fixture
def manager(downloader, storage, mocker):
    tempdir = tempfile.mkdtemp()
    manager = DataManager(Cache(downloader, storage, tempdir))
    manager._tempdir = tempdir
    m = mock.Mock()
    m.headers = {'Content-Disposition': 'test_file'}
    mocker.patch('sunpy.data.data_manager.cache.urlopen', return_value=m)
    yield manager
    shutil.rmtree(tempdir)


@pytest.fixture
def data_function(manager):
    @manager.require('test_file', ['url1/test_file', 'url2'], mocks.MOCK_HASH)
    def foo(manager_tester=lambda x: 1):
        manager_tester(manager)

    return foo
