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
        pass

    def download(self, url, path):
        downloader = Downloader()
        filename = path.split('/')[-1]
        directory = path[:-len(filename)]
        downloader.enqueue_file(url, directory, filename)
        downloader.download()
