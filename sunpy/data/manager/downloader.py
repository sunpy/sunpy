from abc import ABCMeta, abstractmethod


class DownloaderBase(metaclass=ABCMeta):
    @abstractmethod
    def download(self, url, path):
        """
        Downloads and returns the path
        """
        raise NotImplementedError


class MockDownloader(DownloaderBase):
    """MockDownloader"""

    def __init__(self):
        pass

    def download(self, url, path):
        return "/tmp/lol"
