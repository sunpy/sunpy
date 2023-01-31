from unittest.mock import patch

import parfive
import pytest
from packaging.version import Version
from parfive import Results

from sunpy.data.data_manager.downloader import DownloaderError, ParfiveDownloader


def test_ParfiveDownloader_errors():
    """
    Test that ParfiveDownloader raises an error when the download fails.
    """

    downloader = ParfiveDownloader()
    results = Results()

    # TODO: Remove this when we support parfive 2.0.
    if Version(parfive.__version__) >= Version("2.0.0"):
        from parfive.results import Error
    else:
        Error = results._error
    results.errors.append(Error("", "FAKE_URL", ValueError("TEST_ERROR")))
    with patch('parfive.Downloader.download') as download:
        download.return_value = results
        with pytest.raises(DownloaderError, match='TEST_ERROR'):
            downloader.download('https://www.fakewebsite.com', 'test_file')
