from sunpy.data._sample import download_sample_data
from sunpy.data.data_manager.cache import Cache
from sunpy.data.data_manager.storage import SqliteStorage
from sunpy.data.data_manager.downloader import ParfiveDownloader
from sunpy.data.data_manager.manager import DataManager
from sunpy.util.config import get_and_create_download_dir

_download_dir = get_and_create_download_dir()
manager = DataManager(
    Cache(
        ParfiveDownloader(),
        SqliteStorage(_download_dir + '/data_manager.db'),
        _download_dir
    )
)
__all__ = ["download_sample_data", "manager"]
