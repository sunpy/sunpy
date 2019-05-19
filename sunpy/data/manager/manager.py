from contextlib import contextmanager


class DataManager:
    """
    DataManager
    """

    def __init__(self, downloader, storage):
        self._downloader = downloader
        self._storage = storage

        self._skip_hash_check = False
        self._skip_file_names = []

    def require(self, func, name, urls, sha_hash):
        def function(*args, **kwargs):
            # skip hash check = down
            # skip file name = down or replace
            # check store
            if name in self._skip_file_names:
                file_path, _ = self._download_and_hash(urls)
            else:
                details = self._storage.find_by_hash(sha_hash)
                if not details:
                   file_path, file_hash =


            return func(*args, **kwargs)

        return function

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

    def _download_and_hash(self, urls):
        # TODO: Handle multiple urls
        path = self._downloader.download(urls[0])

        # TODO: Calculate the hash
        shahash = 'asdf'

        return path, shahash
