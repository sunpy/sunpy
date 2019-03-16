import io
import os
import sys

from sunpy import config
from sunpy.util.config import (get_and_create_sample_dir, get_and_create_download_dir,
                               CONFIG_DIR, print_config,
                               _find_config_files, _get_user_configdir, _is_writable_dir)

USER = os.path.expanduser('~')


def test_is_writable_dir(tmpdir, tmp_path):
    writeable_dir = tmpdir.mkdir("can_you_right_to_me")
    just_a_path = tmp_path / "sub"
    assert _is_writable_dir(writeable_dir)
    assert not _is_writable_dir(just_a_path)


def test_get_user_configdir(tmpdir):
    # Default
    assert USER in CONFIG_DIR
    assert CONFIG_DIR.split(os.path.sep)[-1] == "sunpy"
    assert CONFIG_DIR == _get_user_configdir()
    # Try to set a manual one.
    tmp_config_dir = tmpdir.mkdir("config_here")
    os.environ["SUNPY_CONFIGDIR"] = tmp_config_dir.strpath
    assert tmp_config_dir == _get_user_configdir()
    # Bypass this under windows.
    if not (os.name == "nt"):
        os.unsetenv("SUNPY_CONFIGDIR")
    del os.environ["SUNPY_CONFIGDIR"]


def test_print_config_files():
    # TODO: Tidy this up.
    stdout = sys.stdout
    out = io.StringIO()
    sys.stdout = out
    print_config()
    sys.stdout = stdout
    out.seek(0)
    printed = out.read()
    assert "time_format = %Y-%m-%d %H:%M:%S" in printed
    assert _find_config_files()[0] in printed
    assert get_and_create_download_dir() in printed
    assert get_and_create_sample_dir() in printed


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
