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

        Parameters
        ----------
        urls : `list` of path-like or one path-like
            A list of urls or a single url.
            The list is for urls of the same file but from different sources.
        namespace : `str`, optional
            A namespace to be used for the file name.
            Defaults to an empty string.
        redownload : `bool`, optional
            Whether to skip cache and redownload.
            Defaults to `False`.

        Returns
        -------
        `pathlib.Path`
            Path to the downloaded file.
        """
        if isinstance(urls, str | Path):
            urls = [urls]
        # Logic plan
        # 1. Check if the file is present in cache by url
        # 2. If present and it has not expired nor redownload, return the file path
        # 3. If not present or present and (expired or redownload), download the file and update cache
        #   a. If there is an error from the above steps, we will return the file from the cache if present
        for url in urls:
            cache_details = self._get_by_url(url)
            if cache_details:
                break
        # If we have a cache hit and we do not want to redownload it and its still valid
        # We want to just return the file
        if cache_details and not redownload and not self._has_expired(cache_details):
            return Path(cache_details['file_path'])
        try:
            file_path, file_hash, url = self._download_and_hash(urls, namespace)
            # We explicitly check for redownload or expiry:
            # This isn't needed as it should be caught by the above if statement
            # but this is more explicit
            if cache_details and (redownload or self._has_expired(cache_details)):
                self._storage.delete_by_key('url', cache_details['url'])
            self._storage.store({
                'file_hash': file_hash,
                'file_path': str(file_path),
                'url': url,
                'time': datetime.now().isoformat(),
             })
            return file_path
        except Exception as e:
            if not cache_details:
                raise e
            exception_msg = f"{e!r} \n"
            warn_user(
                f"{exception_msg}Due to the above error, you will be working with a stale version of the file in the cache."
            )
            return Path(cache_details['file_path'])

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
        return self._storage.find_by_key('file_hash', sha_hash)

    def _get_by_url(self, url):
        """
        Returns the details which is matched by url if present in cache.

        Parameters
        ----------
        url : `str`
            URL of the file.
        """
        return self._storage.find_by_key('url', url)

    def _download_and_hash(self, urls, namespace):
        """
        Downloads the file and returns the path, hash, and URL it used to download.

        Parameters
        ----------
        urls : `list` of `str`
            List of URLs.

        Returns
        -------
        `pathlib.Path`, `str`, `str`
            Path to the downloaded file, SHA-256 hash, and the URL used.
        """
        errors = []
        for url in urls:
            try:
                path = self._cache_dir / (namespace + get_filename(urlopen(url), url))

                self._downloader.download(url, path, overwrite=True)
                shahash = hash_file(path)
                return path, shahash, url
            except Exception as e:
                errors.append(e)
        else:
            if len(errors) == 1:
                raise errors[0]
            msg = "Download failed for all URLs, the first error is shown above."
            raise RuntimeError(msg) from errors[0]
