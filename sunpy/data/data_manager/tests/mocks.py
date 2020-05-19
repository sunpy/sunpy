from sunpy.data.data_manager.downloader import DownloaderBase

# This is the hash of file containing just the character 'a'
MOCK_HASH = "ca978112ca1bbdcafac231b39a23dc4da786eff8147c4e72b9807785afee48bb"


class MockDownloader(DownloaderBase):
    """
    MockDownloader.
    """

    def __init__(self):
        self.times_called = 0
        self.last_called_url = ''

    def download(self, url, path):
        write_to_test_file(path, "a")
        self.times_called += 1
        self.last_called_url = url
        return path


def write_to_test_file(path, contents):
    with open(path, 'w') as f:
        f.write(contents)
