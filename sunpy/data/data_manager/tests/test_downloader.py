from unittest.mock import patch

import pytest
from parfive import Results
from parfive.results import Error

from sunpy.data.data_manager.downloader import DownloaderError, ParfiveDownloader


def test_ParfiveDownloader_errors():
    """
    Test that ParfiveDownloader raises an error when the download fails.
    """
    downloader = ParfiveDownloader()
    results = Results()
    results.errors.append(Error("", "FAKE_URL", ValueError("TEST_ERROR")))
    with patch('parfive.Downloader.download') as download:
        download.return_value = results
        with pytest.raises(DownloaderError, match='TEST_ERROR'):
            downloader.download('https://www.fakewebsite.com', 'test_file')
