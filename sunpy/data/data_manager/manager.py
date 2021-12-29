import pathlib
import functools
from contextlib import contextmanager
from urllib.parse import urlparse

from sunpy.util.exceptions import warn_user
from sunpy.util.util import hash_file

__all__ = ['DataManager']


class DataManager:
    """
    This class provides a remote data manager for managing remote files.

    Parameters
    ----------
    cache: `sunpy.data.data_manager.cache.Cache`
        Cache object to be used by `~sunpy.data.data_manager.manager.DataManager`.
    """

    def __init__(self, cache):
        self._cache = cache

        self._file_cache = {}

        self._namespace = None
        self._skip_hash_check = False
        self._skip_file = {}

    def require(self, name, urls, sha_hash):
        """
        Decorator for informing the data manager about the requirement of
        a file by a function.

        Parameters
        ----------
        name : `str`
            The name to reference the file with.
        urls : `list` or `str`
            A list of urls to download the file from.
        sha_hash : `str`
            SHA-256 hash of file.
        """
        if isinstance(urls, str):
            urls = [urls]

        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                self._namespace = self._get_module(func)
                replace = self._skip_file.get(name)
                if replace:
                    uri_parse = urlparse(replace['uri'])
                    if uri_parse.scheme in ("", "file"):
                        # If a relative file uri is specified (i.e.
                        # `file://sunpy/test`) this maintains compatibility
                        # with the original behaviour where this would be
                        # interpreted as `./sunpy/test` if no scheme is
                        # specified netloc will be '' by default.
                        file_path = uri_parse.netloc + uri_parse.path
                        file_hash = hash_file(file_path)
                    else:
                        file_path, file_hash, _ = self._cache._download_and_hash(
                            [replace['uri']], self._namespace
                        )
                    if replace['hash'] and file_hash != replace['hash']:
                        # if hash provided to replace function doesn't match the hash of the file
                        # raise error
                        raise ValueError(
                            "Hash provided to override_file does not match hash of the file.")
                elif self._skip_hash_check:
                    file_path = self._cache.download(urls, self._namespace, redownload=True)
                else:
                    details = self._cache.get_by_hash(sha_hash)
                    if not details:
                        # In case we are matching by hash and file does not exist
                        # That might mean the wrong hash is supplied to decorator
                        # We match by urls to make sure that is not the case
                        if self._cache_has_file(urls):
                            # If we can't find a file matching sha_hash, but the url is already
                            # in the database
                            raise ValueError(f"{urls} has already been downloaded, but no file "
                                             f"matching the hash {sha_hash} can be found.")
                        file_path = self._cache.download(urls, self._namespace)
                        file_hash = hash_file(file_path)
                        if file_hash != sha_hash:
                            # the hash of the file downloaded does not match provided hash
                            # this means the file has changed on the server.
                            # the function should be updated to use the new
                            # hash. Raise an error to notify.
                            raise RuntimeError(
                                f"Hash of local file ({file_hash}) does not match expected hash ({sha_hash}). "
                                "File may have changed on the remote server.")
                    else:
                        # This is to handle the case when the local file
                        # appears to be tampered/corrupted
                        if hash_file(details['file_path']) != details['file_hash']:
                            warn_user("Hashes do not match, the file will be redownloaded "
                                      "(could be be tampered/corrupted)")
                            file_path = self._cache.download(urls, self._namespace, redownload=True)
                            # Recheck the hash again, if this fails, we will exit.
                            if hash_file(file_path) != details['file_hash']:
                                raise RuntimeError("Redownloaded file also has the incorrect hash."
                                                   "The remote file on the server might have changed.")
                        else:
                            file_path = details['file_path']

                if name not in self._file_cache:
                    self._file_cache[name] = {}
                self._file_cache[name][self._namespace] = file_path
                result = func(*args, **kwargs)
                self._namespace = None
                return result
            return wrapper

        return decorator

    @contextmanager
    def override_file(self, name, uri, sha_hash=None):
        """
        Replaces the file by the name with the file provided by the url/path.

        Parameters
        ----------
        name : `str`
            Name of the file provided in the `require` decorator.
        uri : `str`
            URI of the file which replaces original file. Scheme should be one
            of ``http``, ``https``, ``ftp`` or ``file``. If no scheme is given
            the uri will be interpreted as a local path. i.e.
            ``file:///tmp/test`` and ``/tmp/test`` are the same.
        sha_hash : `str`, optional
            SHA256 hash of the file to compared to after downloading.
        """
        try:
            self._skip_file[name] = {
                'uri': uri,
                'hash': sha_hash,
            }
            yield
        finally:
            _ = self._skip_file.pop(name, None)

    @contextmanager
    def skip_hash_check(self):
        """
        Disables hash checking temporarily

        Examples
        --------
        >>> with remote_data_manager.skip_hash_check():  # doctest: +SKIP
        ...     myfunction()  # doctest: +SKIP
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
        name : `str`
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
        return pathlib.Path(self._file_cache[name][self._namespace])

    def _cache_has_file(self, urls):
        for url in urls:
            if self._cache._get_by_url(url):
                return True
        return False

    def _get_module(self, func):
        """
        Returns the name of the module (appended with a dot) that a function belongs to.

        Parameters
        ----------
        func : function
            A function whose module is to be found.

        Returns
        -------
        A `str` that represents the module name appended with a dot.
        """
        module = func.__module__.lstrip('.').split('.')[0] + '.'
        if module == '__main__.':
            module = ''
        return module
