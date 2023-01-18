from abc import ABCMeta, abstractmethod
from pathlib import Path

from sunpy.util.parfive_helpers import Downloader

__all__ = ['DownloaderBase', 'DownloaderError', 'ParfiveDownloader']


class DownloaderBase(metaclass=ABCMeta):
    """
    Base class for remote data manager downloaders.
    """
    @abstractmethod
    def download(self, url, path):
        """
        Downloads a file.

        Parameters
        ----------
        url : `str`
            URL of the file to be downloaded.
        path : `pathlib.Path` or `str`
            Path where the file should be downloaded to.

        Raises
        ------
        `DownloaderError`
            DownloaderError is raised when download errors.
        """


class DownloaderError(Exception):
    """
    Error to be raised when a download fails.
    """


class ParfiveDownloader(DownloaderBase):
    """
    Concrete implementation of `~sunpy.data.data_manager.downloader.DownloaderBase`
    using :mod:`parfive`.
    """

    def download(self, url, path):
        downloader = Downloader()
        path = Path(path)
        filename = path.name
        directory = path.parent
        downloader.enqueue_file(url, directory, filename)
        try:
            output = downloader.download()
        except Exception as e:
            raise DownloaderError from e
        if output.errors:
            raise DownloaderError(output.errors[0].exception)
