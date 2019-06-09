from abc import ABCMeta, abstractmethod
from parfive import Downloader


class DownloaderBase(metaclass=ABCMeta):
    @abstractmethod
    def download(self, url, path):
        """
        Downloads and returns the path
        """
        raise NotImplementedError


class ParfiveDownloader(DownloaderBase):
    def __init__(self):
        self.downloader = Downloader()

    def download(self, url, path):
        self.downloader.enqueue_file(url, path)
        self.downloader.download()
