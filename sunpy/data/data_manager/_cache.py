import os
from pathlib import Path
from datetime import datetime
from urllib.request import urlopen

import astropy.units as u
from astropy.time import TimeDelta

from sunpy.util.exceptions import warn_user
from sunpy.util.net import get_filename
from sunpy.util.util import hash_file

__all__ = ['Cache']


class Cache:
    """
    Cache provides a way to download and cache files.

    Parameters
    ----------
    downloader: Implementation of `~sunpy.data.data_manager.downloader.DownloaderBase`
        Downloader object for downloading remote files.
    storage: Implementation of `~sunpy.data.data_manager.storage.StorageProviderBase`
        Storage to store metadata about the files.
    cache_dir: `str` or `pathlib.Path`
        Directory where the downloaded files will be stored.
    expiry: `astropy.units.quantity.Quantity` or `None`, optional
        The interval after which the cache is invalidated. If the expiry is `None`,
        then the expiry is not checked. Defaults to 10 days.
    """

    def __init__(self, downloader, storage, cache_dir, expiry=10*u.day):
        self._downloader = downloader
        self._storage = storage
        self._cache_dir = Path(cache_dir)
        self._expiry = expiry if expiry is None else TimeDelta(expiry)

    def download(self, urls, namespace='', redownload=False):
        """
        Downloads the files from the urls.

        The overall flow of this function is:
            1. If ``redownload``: Download, update cache and return file path.
            2. If not ``redownload``: Check cache,
                i. If present in cache:
                    - If cache has expired, remove the entry from cache, download and add to cache
                    - If cache has not expired, return the path

        Parameters
        ----------
        urls : `list` or `str`
            A list of urls or a single url.
        redownload : `bool`
            Whether to skip cache and redownload.

        Returns
        -------
        `pathlib.PosixPath`
            Path to the downloaded file.
        """
        if isinstance(urls, str):
            urls = [urls]
        # Program flow
        # 1. If redownload: Download, update cache and return file path
        # 2. If not redownload: Check cache,
        #    i. If present in cache:
        #        - If cache expired, remove entry from cache, download and add to cache
        #        - If cache not expired, return path
        details = None
        for url in urls:
            details = self._get_by_url(url)
            if details:
                break
        if details:
            if redownload or self._has_expired(details):
                # if file is in cache and it has to be redownloaded or the cache has expired
                # then remove the file and delete the details from the storage
                os.remove(details['file_path'])
                self._storage.delete_by_key('url', details['url'])
            else:
                return Path(details['file_path'])

        file_path, file_hash, url = self._download_and_hash(urls, namespace)

        self._storage.store({
            'file_hash': file_hash,
            'file_path': str(file_path),
            'url': url,
            'time': datetime.now().isoformat(),
        })
        return file_path

    def _has_expired(self, details):
        """
        Whether the url corresponding to details in cache has expired or not.

        Parameters
        ----------
        details : `dict`
            Details detached from cache.

        Returns
        -------
        `bool`
            Whether the url has expired or not.
        """
        time = details.get("time", datetime.now().isoformat())
        time = datetime.fromisoformat(time)
        return self._expiry and datetime.now() - time > self._expiry

    def get_by_hash(self, sha_hash):
        """
        Returns the details which is matched by hash if present in cache.

        Parameters
        ----------
        sha_hash : `str`
            SHA-256 hash of the file.
        """
        details = self._storage.find_by_key('file_hash', sha_hash)
        return details

    def _get_by_url(self, url):
        """
        Returns the details which is matched by url if present in cache.

        Parameters
        ----------
        url : `str`
            URL of the file.
        """
        details = self._storage.find_by_key('url', url)
        return details

    def _download_and_hash(self, urls, namespace=''):
        """
        Downloads the file and returns the path, hash and url it used to download.

        Parameters
        ----------
        urls : `list`
            List of urls.

        Returns
        -------
        `str`, `str`, `str`
            Path, hash and URL of the file.
        """
        def download(url):
            path = self._cache_dir / (namespace + get_filename(urlopen(url), url))
            self._downloader.download(url, path)
            shahash = hash_file(path)
            return path, shahash, url

        errors = []
        for url in urls:
            try:
                return download(url)
            except Exception as e:
                warn_user(f"{e}")
                errors.append(f"{e}")
        else:
            raise RuntimeError(errors)
