import os

from sunpy import config
from sunpy.util.config import get_and_create_sample_dir, get_and_create_download_dir

USER = os.path.expanduser('~')


def test_get_and_create_download_dir():
    # test default config
    path = get_and_create_download_dir()
    assert path == os.path.join(USER, 'sunpy', 'data')
    # test updated config
    new_path = os.path.join(USER, 'sunpy_data_here_please')
    config.set('downloads', 'download_dir', new_path)
    path = get_and_create_download_dir()
    assert path == os.path.join(USER, new_path)
    # Set the config back
    os.rmdir(new_path)
    config.set('downloads', 'download_dir', os.path.join(USER, 'sunpy', 'data'))


def test_get_and_create_sample_dir():
    # test default config
    path = get_and_create_sample_dir()
    assert path == os.path.join(USER, 'sunpy', 'data', 'sample_data')
    # test updated config
    new_path = os.path.join(USER, 'sample_data_here_please')
    config.set('downloads', 'sample_dir', new_path)
    path = get_and_create_sample_dir()
    assert path == new_path
    # Set the config back
    os.rmdir(new_path)
    config.set('downloads', 'sample_dir', os.path.join(USER, 'sunpy', 'data', 'sample_data'))
