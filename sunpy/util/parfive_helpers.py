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

        # Only specify headers to parfive 1.1
        if not parfive.__version__.startswith('1.0'):
            kwargs["headers"] = headers

        super().__init__(*args, **kwargs)
