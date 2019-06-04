class Cache:
    def __init__(self, downloader, storage):
        self._downloader = downloader
        self._storage = storage

    def download(self, urls, redownload=False):
        # TODO: Expiry time
        # XXX: Expiry time cache level or download level?
        if not redownload:
            details = self._get_by_url(urls[0])
            if details:
                return details['file_path']

        file_path, file_hash, url = self._download_and_hash(urls)

        if not redownload:
            self._storage.store({
                'file_hash': file_hash,
                'file_path': file_path,
                'url': url,
            })
        return file_path

    def get_by_hash(self, sha_hash):
        details = self._storage.find_by_key('file_hash', sha_hash)
        return details

    def _get_by_url(self, url):
        details = self._storage.find_by_key('url', url)
        return details

    def _download_and_hash(self, urls):
        # TODO: Handle multiple urls
        path = self._downloader.download(urls[0])

        # TODO: Calculate the hash
        shahash = 'asdf'

        return path, shahash, urls[0]
