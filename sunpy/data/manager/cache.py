from urllib.request import urlopen

from sunpy.util.util import hash_file, replacement_filename
from sunpy.util.net import get_filename


class Cache:
    """
    Cache handles caching.
    """

    def __init__(self, downloader, storage, cache_dir):
        self._downloader = downloader
        self._storage = storage
        self._cache_dir = cache_dir
        if not self._cache_dir.endswith('/'):
            self._cache_dir += '/'

    def download(self, urls, redownload=False):
        """
        Downloads the files from the urls.

        Parameters
        ----------
        urls: `list`
            A list of urls.
        redownload: `bool`
            Whether to skip cache and redownload.
        """
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
        """
        Returns the details which is matched by hash if present in cache.

        Parameters
        ----------
        sha_hash: `str`
            SHA-1 hash of the file.
        """
        details = self._storage.find_by_key('file_hash', sha_hash)
        return details

    def _get_by_url(self, url):
        """
        Returns the details which is matched by url if present in cache.

        Parameters
        ----------
        url: `str`
            URL of the file.
        """
        details = self._storage.find_by_key('url', url)
        return details

    def _download_and_hash(self, urls):
        """
        Downloads the file and returns the path, hash and url it used to download.

        Parameters
        ----------
        urls: `list`
            List of urls.

        Returns
        -------
        `str`, `str`, `str`
            Path, hash and URL of the file.
        """
        # TODO: Handle multiple urls
        url = urls[0]
        path = self._cache_dir + get_filename(urlopen(url), url)
        path = replacement_filename(path)
        self._downloader.download(url, path)

        shahash = hash_file(path)

        return path, shahash, urls[0]
