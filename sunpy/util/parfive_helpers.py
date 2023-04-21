import os
import sys
from functools import wraps

import aiohttp
import parfive
from packaging.version import Version
from parfive import Results, SessionConfig

import sunpy

__all__ = ["Downloader", "Results"]

parfive_version = Version(parfive.__version__)
sunpy_headers = {
    "User-Agent": f"sunpy/{sunpy.__version__} parfive/{parfive.__version__} "
    f"aiohttp/{aiohttp.__version__} python/{sys.version[:5]} "
}
if os.environ.get("SUNPY_PYTEST_RUN"):
    # We want to allow remote servers to filter users from our online tests
    sunpy_headers["User-Agent"] = f"sunpy-CI-tests/{sunpy.__version__} parfive/{parfive.__version__} "
config = SessionConfig(headers=sunpy_headers)


class Downloader(parfive.Downloader):
    @wraps(parfive.Downloader.__init__)
    def __init__(self, *args, **kwargs):
        if "config" not in kwargs:
            kwargs["config"] = config
        super().__init__(*args, **kwargs)
