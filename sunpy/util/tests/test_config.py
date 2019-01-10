from pathlib import Path

from sunpy import config
from sunpy.util.config import get_and_create_sample_dir, get_and_create_download_dir

USER = Path.home()


def test_get_and_create_download_dir():
    # test default config
    path = get_and_create_download_dir()
    assert path == str(Path.home().joinpath(USER, 'sunpy', 'data'))
    # test updated config
    new_path = str(Path.home().joinpath(USER, 'data_here_please'))
    config.set('downloads', 'download_dir', new_path)
    path = get_and_create_download_dir()
    assert path == str(Path.home().joinpath(USER, new_path))


def test_get_and_create_sample_dir():
    # test default config
    path = get_and_create_sample_dir()
    assert path == str(Path.home().joinpath(USER, 'sunpy', 'data', 'sample_data'))
    # test updated config
    new_path = str(Path.home().joinpath(USER, 'data_here_please'))
    config.set('downloads', 'sample_dir', new_path)
    path = get_and_create_sample_dir()
    assert path == new_path
