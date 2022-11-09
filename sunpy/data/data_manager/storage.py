"""
Storage module contains the abstract implementation of storage
for `sunpy.data.data_manager.Cache` and a concrete implementation
using sqlite.
"""
import sqlite3
from abc import ABCMeta, abstractmethod
from pathlib import Path
from contextlib import contextmanager

__all__ = [
    'StorageProviderBase',
    'SqliteStorage',
    'InMemStorage',
]


class StorageProviderBase(metaclass=ABCMeta):
    """
    Base class for remote data manager storage providers.
    """
    @abstractmethod
    def find_by_key(self, key, value):
        """
        Returns the file details if value corresponding to the key
        found in storage. Returns `None` if hash not found.

        Parameters
        ----------
        key : `str`
            The key/column name of the field.
        value : `str`
            The value associated with the key of the entry.

        Returns
        -------
        `dict` or `None`
            `dict` contains the details of the file. `None` if hash not found.

        Raises
        ------
        ``KeyError``
            KeyError is raised if key does not exist.
        """

    @abstractmethod
    def delete_by_key(self, key, value):
        """
        Deletes the matching entry from the store.

        Parameters
        ----------
        key : `str`
            The key/column name of the field.
        value : `str`
            The value associated with the key of the entry.

        Raises
        ------
        ``KeyError``
            KeyError is raised if key does not exist.
        """

    @abstractmethod
    def store(self, details):
        """
        Stores the details in the storage.

        Parameters
        ----------
        details : `dict`
            Details to be stored.
        """


class InMemStorage(StorageProviderBase):
    """
    This provides a storage stored in memory.
    """

    def __init__(self):
        self._store = []

    def store(self, details):
        self._store += [details]

    def delete_by_key(self, key, value):
        for i in self._store:
            if i[key] == value:
                self._store.remove(i)

    def find_by_key(self, key, value):
        for i in self._store:
            if i[key] == value:
                return i
        return None


class SqliteStorage(StorageProviderBase):
    """
    This provides a sqlite backend for storage.

    Parameters
    ----------
    path: `str`
        Path to the database file.
    """
    COLUMN_NAMES = [
        'file_hash',
        'file_path',
        'url',
        'time',
    ]

    def __init__(self, path):
        self._db_path = Path(path)
        self._table_name = 'cache_storage'

        self._db_path.parent.mkdir(parents=True, exist_ok=True)
        if not self._db_path.exists():
            # setup database
            self._setup()

    def _setup(self):
        with self.connection(commit=True) as conn:
            self._create_table(conn)

    def _create_table(self, conn):
        schema = ' text, '.join(self.COLUMN_NAMES) + ' text'
        conn.execute(f'''CREATE TABLE IF NOT EXISTS {self._table_name} ({schema})''')

    @contextmanager
    def connection(self, commit=False):
        """
        A context manager which provides an easy way to handle db connections.

        Parameters
        ----------
        commit : `bool`
            Whether to commit after successful execution of db command.
        """
        conn = sqlite3.connect(str(self._db_path))
        self._create_table(conn)
        try:
            yield conn
            if commit:
                conn.commit()
        finally:
            conn.close()

    def find_by_key(self, key, value):
        if key not in self.COLUMN_NAMES:
            raise KeyError
        with self.connection() as conn:
            cursor = conn.cursor()
            cursor.execute(f'''SELECT * FROM {self._table_name}
                                      WHERE {key}="{value}"''')
            row = cursor.fetchone()
            if row:
                return dict(zip(self.COLUMN_NAMES, row))
            return None

    def delete_by_key(self, key, value):
        if key not in self.COLUMN_NAMES:
            raise KeyError
        with self.connection(commit=True) as conn:
            cursor = conn.cursor()
            cursor.execute(f'''DELETE FROM {self._table_name}
                                      WHERE {key}="{value}"''')

    def store(self, details):
        values = [details[k] for k in self.COLUMN_NAMES]
        placeholder = '?,' * len(values)
        placeholder = placeholder[:-1]
        with self.connection(commit=True) as conn:
            conn.execute(f'''INSERT INTO {self._table_name}
                             VALUES ({placeholder})''', list(values))
