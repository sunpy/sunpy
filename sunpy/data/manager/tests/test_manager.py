from sunpy.data.manager.manager import DataManager
from sunpy.data.manager.storage import InMemStorage
from sunpy.data.manager.downloader import DownloaderBase


class MockDownloader(DownloaderBase):
    """MockDownloader"""

    def __init__(self):
        self.times_called = 0

    def download(self, url):
        self.times_called += 1
        return "/tmp/lol"


def test_basic():
    downloader = MockDownloader()
    storage = InMemStorage()
    manager = DataManager(downloader, storage)

    @manager.require('test_file', ['url1', 'url2'], 'hash')
    def foo():
        pass
    foo()

    assert downloader.times_called == 1
