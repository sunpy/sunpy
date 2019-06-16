from abc import ABCMeta, abstractmethod
from parfive import Downloader
from pathlib import Path


class DownloaderBase(metaclass=ABCMeta):
    @abstractmethod
    def download(self, url, path):
        """
        Downloads and returns the path
        """
        raise NotImplementedError


class ParfiveDownloader(DownloaderBase):
    def __init__(self):
        pass

    def download(self, url, path):
        downloader = Downloader()
        path = Path(path)
        filename = path.name
        directory = path.parent
        downloader.enqueue_file(url, directory, filename)
        downloader.download()
