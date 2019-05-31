from contextlib import contextmanager

import functools


class DataManager:
    """
    DataManager
    """

    def __init__(self, downloader, storage):
        self._downloader = downloader
        self._storage = storage

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
                # skip hash check = down
                # skip file name = down or replace
                # check store
                replace = self._skip_file.get(name, None)
                if replace:
                    pass
                elif self._skip_hash_check:
                    file_path, _ = self._download_and_hash(urls)
                else:
                    details = self._storage.find_by_hash(sha_hash)
                    if not details:
                        file_path, file_hash = self._download_and_hash(urls)
                        # TODO: add modified date
                        # TODO: add function names
                        self._storage.store({
                            'file_hash': file_hash,
                            'file_path': file_path,
                        })
                    else:
                        file_path, file_hash = details['file_path'], details['file_hash']

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
        file: str
        Path of the file

        Raises
        ------
        KeyError
        If name is not in the cache

        TODO: Return `Path` instead?
        """
        return self._file_cache[name]

    def _download_and_hash(self, urls):
        # TODO: Handle multiple urls
        # TODO: Calculate path here. Don't depened on subclasses of downloader to use unique path
        path = self._downloader.download(urls[0])

        # TODO: Calculate the hash
        shahash = 'asdf'

        return path, shahash
