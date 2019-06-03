class Cache:
    def __init__(self, downloader, storage):
        self._downloader = downloader
        self._storage = storage

    def download(self, urls, redownload=False, return_hash=False):
        # TODO: Cache
        file_path, file_hash, url = self._download_and_hash(urls)
        self._storage.store({
            'file_hash': file_hash,
            'file_path': file_path,
            'url': url,
        })
        if return_hash:
            return file_path, file_hash
        return file_path

    def get_by_hash(self, sha_hash):
        details = self._storage.find_by_key('file_hash', sha_hash)
        return details['file_path']

    def _get_by_url(self, url):
        details = self._storage.find_by_key('url', url)
        return details['file_path']

    def _download_and_hash(self, urls):
        # TODO: Handle multiple urls
        path = self._downloader.download(urls[0])

        # TODO: Calculate the hash
        shahash = 'asdf'

        return path, shahash, urls[0]
