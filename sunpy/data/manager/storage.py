from abc import ABCMeta, abstractmethod


class StorageProviderBase(metaclass=ABCMeta):
    """Base class for remote data manager storage providers"""
    @abstractmethod
    def find_by_key(self, key, value):
        """
        Returns the file details if hash found in storage.
        Returns `None` if hash not found.

        Parameters
        ----------
        file_hash: `str`
        Hash of the file

        Returns
        -------
        `dict` or `None`
        `dict` contains the details of the file. `None` if hash not found.
        """
        raise NotImplementedError

    @abstractmethod
    def store(self, details):
        """
        Stores the details in the storage.

        Parameters
        ----------
        details: `dict`
        Details to be stored
        """
        raise NotImplementedError


class InMemStorage(StorageProviderBase):
    def __init__(self):
        self._store = []

    def store(self, details):
        self._store += [details]

    def find_by_key(self, key, value):
        for i in self._store:
            if i['file_hash'] == value:
                return i
        return None
