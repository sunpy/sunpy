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
        self._require_files = {}

    def require(self, name, urls, sha_hash, defer_download=False):
        """
        Decorator for informing the data manager about the requirement of
        a file by a function. Optionally defer downloading the file until it's
        requested in `~sunpy.data.data_manager.manager.DataManager.get`.

        Parameters
        ----------
        name : `str`
            The name to reference the file with.
        urls : `list` or `str`
            A list of urls to download the file from.
        sha_hash : `str`
            SHA-256 hash of file.
        defer_download : `bool`, optional
            If `True`, the file download is deferred until it is requested via `~sunpy.data.data_manager.manager.DataManager.get`.
        """
        if isinstance(urls, str):
            urls = [urls]

        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                self._namespace = self._get_module(func)
                if defer_download:
                    self._require_files[name] = {
                        'urls': urls,
                        'sha_hash': sha_hash,
                    }
                else:
                    self._download_and_cache_file(name, urls, sha_hash)

                result = func(*args, **kwargs)
                self._namespace = None
                return result
            return wrapper

        return decorator

    def _download_and_cache_file(self, name, urls, sha_hash):
        """
        Internal method to handle the downloading and caching logic.
        """
        replace = self._skip_file.get(name)
        if replace:
            uri_parse = urlparse(replace['uri'])
            if uri_parse.scheme in ("", "file"):
                file_path = uri_parse.netloc + uri_parse.path
                file_hash = hash_file(file_path)
            else:
                file_path, file_hash, _ = self._cache._download_and_hash([replace['uri']], self._namespace)
            if replace['hash'] and file_hash != replace['hash']:
                raise ValueError("Hash provided to override_file does not match hash of the file.")
        elif self._skip_hash_check:
            file_path = self._cache.download(urls, self._namespace, redownload=True)
        else:
            details = self._cache.get_by_hash(sha_hash)
            if not details:
                if self._cache_has_file(urls):
                    raise ValueError(f"{urls} has already been downloaded, but no file matching the hash {sha_hash} can be found.")
                file_path = self._cache.download(urls, self._namespace)
                file_hash = hash_file(file_path)
                if file_hash != sha_hash:
                    raise RuntimeError(f"Hash of local file ({file_hash}) does not match expected hash ({sha_hash}). File may have changed on the remote server.")
            else:
                if not pathlib.Path(details["file_path"]).is_file():
                    warn_user("Requested file appears to missing and will be redownloaded.")
                    self._cache._download_and_hash(urls, self._namespace)

                if hash_file(details['file_path']) != details['file_hash']:
                    warn_user("Hashes do not match, the file will be redownloaded (could be tampered/corrupted)")
                    file_path = self._cache.download(urls, self._namespace, redownload=True)
                    if hash_file(file_path) != details['file_hash']:
                        raise RuntimeError("Redownloaded file also has the incorrect hash. The remote file on the server might have changed.")
                else:
                    file_path = details['file_path']

        if name not in self._file_cache:
            self._file_cache[name] = {}
        self._file_cache[name][self._namespace] = file_path

    def get(self, name):
        """
        Get the file by name, and download it if deferred.

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
        if name in self._require_files:
            file_info = self._require_files.pop(name)
            self._download_and_cache_file(name, file_info['urls'], file_info['sha_hash'])

        return pathlib.Path(self._file_cache[name][self._namespace])

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
        Disables hash checking temporarily.

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
