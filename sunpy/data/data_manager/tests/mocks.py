from sunpy.data.data_manager.downloader import DownloaderBase

MOCK_HASH = "86f7e437faa5a7fce15d1ddcb9eaeaea377667b8"

class MockDownloader(DownloaderBase):
    """
    MockDownloader.
    """

    def __init__(self):
        self.times_called = 0

    def download(self, url, path):
        write_to_test_file(path, "a")
        self.times_called += 1
        return path


def write_to_test_file(path, contents):
    # TODO: Use tempfile. here and other places.
    with open(path, 'w') as f:
        f.write(contents)
