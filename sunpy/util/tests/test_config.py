from pathlib import Path

from sunpy import config
from sunpy.util.config import get_and_create_sample_dir, get_and_create_download_dir

USER = str(Path.home())


def test_get_and_create_download_dir():
    # test default config
    path = get_and_create_download_dir()
    assert path == str(Path(USER).joinpath('sunpy', 'data'))
    # test updated config
    new_path = str(Path(USER).joinpath('data_here_please'))
    config.set('downloads', 'download_dir', new_path)
    path = get_and_create_download_dir()
    assert path == str(Path(USER).joinpath(new_path))

def test_get_and_create_sample_dir():
    # test default config
    path = get_and_create_sample_dir()
    assert path == str(Path(USER).joinpath('sunpy', 'data', 'sample_data'))
    # test updated config
    new_path = str(Path(USER).joinpath('data_here_please'))
    config.set('downloads', 'sample_dir', new_path)
    path = get_and_create_sample_dir()
    assert path == new_path
    # Set the config back
    config.set('downloads', 'download_dir', os.path.join(USER, 'sunpy', 'data', 'sample_data'))
