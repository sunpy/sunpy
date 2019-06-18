import pathlib
import functools
from contextlib import contextmanager

from sunpy.util.util import hash_file


class DataManager:
    """
    DataManager

    Parameters
    ----------
    cache: `sunpy.data.data_manager.cache.Cache`
    """

    def __init__(self, cache):
        self._cache = cache

        self._file_cache = {}

        self._skip_hash_check = False
        self._skip_file = {}  # Dict[str, str]

    def require(self, name, urls, sha_hash):
        """
        Decorator for informing the data manager about the requirement of
        a file by a function.

        Parameters
        ----------
        name: `str`
            The name to reference the file with.
        urls: `list`
            A list of urls to download the file from.
        sha_hash: `str`
            SHA-1 hash of file.
        """
        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                replace = self._skip_file.get(name, None)
                if replace:
                    if replace.startswith('file://'):
                        file_path = replace[len('file://'):]
                    else:
                        file_path = self._cache.download([replace], redownload=True)
                elif self._skip_hash_check:
                    file_path = self._cache.download(urls, redownload=True)
                else:
                    details = self._cache.get_by_hash(sha_hash)
                    err = KeyError("Hash does not match")
                    if not details:
                        # In case we are matching by hash and file does not exist
                        # That might mean the wrong hash is supplied to decorator
                        # We match by urls to make sure that is not the case
                        # XXX: Is this required? Would be useful if users want to
                        # write functions using `data_manager` or if developers
                        # are sleepy while writing the code
                        # TODO: Raise a warning
                        # https://github.com/sunpy/sunpy/pull/3124#discussion_r291576259
                        if self._cache_has_file(urls):
                            raise err
                        file_path = self._cache.download(urls)
                    else:
                        # This is to handle the case when the file is tampered on disk
                        if hash_file(details['file_path']) != details['file_hash']:
                            raise err
                        file_path = details['file_path']

                self._file_cache[name] = file_path
                return func(*args, **kwargs)
            return wrapper

        return decorator

    @contextmanager
    def replace_file(self, name, uri):
        """
        Replaces the file by the name with the file provided by the url/path.

        TODO: Hash

        Parameters
        ----------
        name: `str`
            Name of the file provided in the `require` decorator.
        uri: `str`
            URI of the file which replaces original file. One of ``http``, ``https``, ``ftp``
            or ``file``.
        """
        try:
            self._skip_file[name] = uri
            yield
        finally:
            _ = self._skip_file.pop(name, None)

    @contextmanager
    def skip_hash_check(self):
        """
        Disables hash checking temporarily

        Examples
        --------
        >>> with remote_data_manager.skip_hash_check():
        ...     myfunction()
        """
        try:
            self._skip_hash_check = True
            yield
        finally:
            self._skip_hash_check = False

    def get(self, name):
        """
        Get the file by name.

        Parameters
        ----------
        name: `str`
            Name of the file given to the data manager, same as the one provided
            in `~sunpy.data.data_manager.manager.DataManager.require`.

        Returns
        -------
        `pathlib.Path`
            Path of the file.

        Raises
        ------
        `KeyError`
            If ``name`` is not in the cache.
        """
        return pathlib.Path(self._file_cache[name])

    def _cache_has_file(self, urls):
        for url in urls:
            if self._cache._get_by_url(url):
                return True
        return False
