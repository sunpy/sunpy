import os
import sys
from functools import wraps

import aiohttp
import parfive
from packaging.version import Version
from parfive import Results

import sunpy

__all__ = ["Downloader", "Results"]

parfive_version = Version(parfive.__version__)
sunpy_headers = {
    "User-Agent": f"sunpy/{sunpy.__version__} parfive/{parfive.__version__} "
    f"aiohttp/{aiohttp.__version__} python/{sys.version[:5]}"
}


if parfive_version < Version("2.0a0"):
    # Overload the parfive downloader class to set the User-Agent string
    class Downloader(parfive.Downloader):
        @wraps(parfive.Downloader.__init__)
        def __init__(self, *args, **kwargs):
            headers = kwargs.pop("headers", {})
            kwargs["headers"] = {**sunpy_headers, **headers}

            # Works with conftest to hide the progress bar.
            progress = os.environ.get("PARFIVE_HIDE_PROGESS", None)
            if progress:
                kwargs["progress"] = False

            super().__init__(*args, **kwargs)

else:
    from parfive import SessionConfig

    config = SessionConfig(headers=sunpy_headers)

    class Downloader(parfive.Downloader):
        @wraps(parfive.Downloader.__init__)
        def __init__(self, *args, **kwargs):
            if "config" not in kwargs:
                kwargs["config"] = config

            # This is to make our sample data download code work on 1.x and 2.x
            # when we depend on 2+ we should remove this.
            headers = kwargs.pop("headers", {})
            kwargs["config"].headers = {**kwargs["config"].headers, **headers}

            super().__init__(*args, **kwargs)
