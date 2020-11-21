import os
import sys
from functools import wraps

import aiohttp
import parfive
from parfive import Results

import sunpy

__all__ = ['Downloader', 'Results']


# Overload the parfive downloader class to set the User-Agent string
class Downloader(parfive.Downloader):
    @wraps(parfive.Downloader.__init__)
    def __init__(self, *args, **kwargs):
        headers = kwargs.pop("headers", None)
        if headers is None or 'User-Agent' not in headers:
            headers = {'User-Agent':
                       f"sunpy/{sunpy.__version__} parfive/{parfive.__version__} "
                       f"aiohttp/{aiohttp.__version__} python/{sys.version[:5]}"}
        kwargs["headers"] = headers
        # Works with conftest to hide the progress bar.
        progress = os.environ.get("HIDE_PARFIVE_PROGESS", None)
        if progress:
            kwargs["progress"] = False
        super().__init__(*args, **kwargs)
