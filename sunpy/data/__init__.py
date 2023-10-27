import astropy.units as u

import sunpy
from sunpy import config
from sunpy.data.data_manager.cache import Cache
from sunpy.data.data_manager.downloader import ParfiveDownloader
from sunpy.data.data_manager.manager import DataManager
from sunpy.data.data_manager.storage import SqliteStorage
from sunpy.util.config import CACHE_DIR

_download_dir = config.get('downloads', 'remote_data_manager_dir')


manager = DataManager(
    Cache(
        ParfiveDownloader(),
        SqliteStorage(_download_dir + '/data_manager.db'),
        _download_dir
    )
)
cache = Cache(
    ParfiveDownloader(),
    SqliteStorage(CACHE_DIR + '/cache.db'),
    CACHE_DIR,
    expiry=int(config.get('downloads', 'cache_expiry')) * u.day
)

__all__ = ["manager", "cache"]
