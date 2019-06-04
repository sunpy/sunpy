from contextlib import contextmanager
import pathlib

import functools


class DataManager:
    """
    DataManager
    """

    def __init__(self, cache):
        # TODO: Folder prefix should be init argument
        self._cache = cache

        self._file_cache = {}

        self._skip_hash_check = False
        self._skip_file = {}  # Dict[str, str]

    def require(self, name, urls, sha_hash):
        """decorator for doing stuff

        Parameters
        ----------
        name: str
        The name to reference the file with
        urls: list
        List of urls to download the file from
        sha_hash: str
        Hash of file
        """
        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                # TODO: Refractor into function(s)
                replace = self._skip_file.get(name, None)
                if replace:
                    if replace.startswith('file://'):
                        file_path = replace[len('file://'):]
                    else:
                        file_path = self._cache.download([replace], redownload=True)
                elif self._skip_hash_check:
                    # XXX: Should this redownload every time
                    file_path = self._cache.download(urls, redownload=True)
                else:
                    details = self._cache.get_by_hash(sha_hash)
                    if not details:
                        file_path = self._cache.download(urls)
                    else:
                        file_path = details['file_path']

                # TODO: handle name clashing
                self._file_cache[name] = file_path
                return func(*args, **kwargs)
            return wrapper

        return decorator

    @contextmanager
    def replace_file(self, name, uri):
        """Replaces the file by the name with the file provided by the url/path

        TODO: Hash

        Parameters
        ----------
        name: str
        Name of the file provided in the `require` decorator
        uri: str
        URI of the file which replaces original file. One of `http`, `https`, `ftp`
        or `file`
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
            with remote_data_manager.skip_hash_check():
                myfunction()
        """
        try:
            self._skip_hash_check = True
            yield
        finally:
            self._skip_hash_check = False

    def get(self, name):
        """get the file by name

        Parameters
        ----------
        name: str
        Name of the file given to the data manager, same as the one provided
        in `~sunpy.data.manager.manager.DataManager.require`

        Returns
        -------
        file: `pathlib.Path`
        Path of the file

        Raises
        ------
        KeyError
        If name is not in the cache
        """
        return pathlib.Path(self._file_cache[name])

    def _download_and_hash(self, urls):
        # TODO: Handle multiple urls
        # TODO: Calculate path here. Don't depened on subclasses of downloader to use unique path
        path = self._downloader.download(urls[0])

        # TODO: Calculate the hash
        shahash = 'hash'

        return path, shahash
