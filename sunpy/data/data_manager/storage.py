import os
import sqlite3
from abc import ABCMeta, abstractmethod
from contextlib import contextmanager


class StorageProviderBase(metaclass=ABCMeta):
    """
    Base class for remote data manager storage providers
    """
    @abstractmethod
    def find_by_key(self, key, value):
        """
        Returns the file details if hash found in storage. Returns `None`
        if hash not found.

        Parameters
        ----------
        file_hash: `str`
            Hash of the file.

        Returns
        -------
        `dict` or `None`
            `dict` contains the details of the file. `None` if hash not found.

        Raises
        ------
        KeyError
             KeyError is raised if key does not exist.
        """
        raise NotImplementedError

    @abstractmethod
    def store(self, details):
        """
        Stores the details in the storage.

        Parameters
        ----------
        details: `dict`
            Details to be stored.
        """
        raise NotImplementedError


class InMemStorage(StorageProviderBase):
    """
    InMemStorage provides a storage stored in memory
    as `dict`s
    """

    def __init__(self):
        self._store = []

    def store(self, details):
        self._store += [details]

    def find_by_key(self, key, value):
        for i in self._store:
            if i[key] == value:
                return i
        return None


class SqliteStorage(StorageProviderBase):
    """
    SqliteStorage provides a sqlite backend for storage.
    """
    COLOUMN_NAMES = [
        'file_hash',
        'file_path',
        'url',
    ]

    def __init__(self, path):
        self._db_path = path
        self._table_name = 'cache_storage'

        if not os.path.exists(self._db_path):
            # setup database
            self._setup()

    def _setup(self):
        schema = ' text, '.join(self.COLOUMN_NAMES) + ' text'
        with self.connection(commit=True) as conn:
            conn.execute(f'''CREATE TABLE {self._table_name}
                             ({schema})''')

    @contextmanager
    def connection(self, commit=False):
        """
        A contextmanger which provides an easy way to handle db connections.

        Parameters
        ----------
        commit: `bool`
        Whether to commit after succesful execution of db command.
        """
        conn = sqlite3.connect(self._db_path)
        try:
            yield conn
            if commit:
                conn.commit()
        finally:
            conn.close()

    def find_by_key(self, key, value):
        if key not in self.COLOUMN_NAMES:
            raise KeyError
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(f'''SELECT * FROM {self._table_name}
                                      WHERE {key}="{value}"''')
            row = cursor.fetchone()
            if row:
                return dict(zip(self.COLOUMN_NAMES, row))
            return None

    def store(self, details):
        values = [details[k] for k in self.COLOUMN_NAMES]
        placeholder = '?,' * len(values)
        placeholder = placeholder[:-1]
        with self.connection(commit=True) as conn:
            conn.execute(f'''INSERT INTO {self._table_name}
                             VALUES ({placeholder})''', list(values))
