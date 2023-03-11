import io
import os
from pathlib import Path
from contextlib import redirect_stdout

import pytest

import sunpy
from sunpy import config
from sunpy.extern.appdirs import AppDirs
from sunpy.util import SunpyUserWarning
from sunpy.util.config import (
    CONFIG_DIR,
    _find_config_files,
    _get_user_configdir,
    _is_writable_dir,
    copy_default_config,
    dirs,
    get_and_create_download_dir,
    get_and_create_sample_dir,
    print_config,
)

USER = os.path.expanduser('~')


def test_is_writable_dir(tmpdir, tmp_path):
    assert _is_writable_dir(tmpdir)
    tmp_file = tmpdir.join("hello.txt")
    # Have to write to the file otherwise its seen as a directory(?!)
    tmp_file.write("content")
    # Checks directory with a file
    assert _is_writable_dir(tmpdir)
    # Checks a filepath instead of directory
    assert not _is_writable_dir(tmp_file)


def test_get_user_configdir(tmpdir, tmp_path, undo_config_dir_patch):
    # Default
    assert USER in CONFIG_DIR
    assert CONFIG_DIR.split(os.path.sep)[-1] == "sunpy"
    assert CONFIG_DIR == _get_user_configdir()
    # Try to set a manual one (already created)
    tmp_config_dir = tmpdir.mkdir("config_here")
    os.environ["SUNPY_CONFIGDIR"] = tmp_config_dir.strpath
    assert tmp_config_dir == _get_user_configdir()
    # Bypass this under windows.
    if not (os.name == "nt"):
        os.unsetenv("SUNPY_CONFIGDIR")
    del os.environ["SUNPY_CONFIGDIR"]
    # Try to set a manual one (not created)
    os.environ["SUNPY_CONFIGDIR"] = str(tmp_path)
    assert str(tmp_path) == _get_user_configdir()
    # Bypass this under windows.
    if not (os.name == "nt"):
        os.unsetenv("SUNPY_CONFIGDIR")
    del os.environ["SUNPY_CONFIGDIR"]


def test_print_config_files(tmpdir, tmp_path, undo_download_dir_patch):
    with io.StringIO() as buf, redirect_stdout(buf):
        print_config()
        printed = buf.getvalue()
    assert "time_format = %Y-%m-%d %H:%M:%S" in printed
    assert _find_config_files()[0] in printed
    assert get_and_create_download_dir() in printed
    assert get_and_create_sample_dir() in printed


def test_find_config_user_site_files(tmpdir, tmp_path, undo_download_dir_patch, monkeypatch):

    # Set a site config dir
    tmp_config_site_dir = tmpdir.mkdir("sunpy_test_configdir_site")
    monkeypatch.setattr(AppDirs, "site_config_dir", tmp_config_site_dir.strpath)

    # Create sunpyrc to site configdir
    tmp_site_file = tmp_config_site_dir.join("sunpyrc")
    # Have to write to the file otherwise its seen as a directory(?!)
    tmp_site_file.write("content")

    # Test that the second config file of two is the site file
    config_files = _find_config_files()
    assert len(config_files) == 2
    assert config_files[1] == tmp_site_file.strpath

    # set a user config dir
    tmp_config_user_dir = tmpdir.mkdir("sunpy_test_configdir")
    monkeypatch.setenv("SUNPY_CONFIGDIR", tmp_config_user_dir.strpath)

    # Copy default sunpyrc to user configdir
    # With a sunpyrc file in the site config dir,
    # the function should raise a warning.
    with pytest.warns(SunpyUserWarning):
        copy_default_config()

    # Test that the second config file of three is the site file
    # and the third of three is the user config file
    config_files = _find_config_files()
    assert len(config_files) == 3
    assert config_files[1] == tmp_site_file.strpath
    assert config_files[2] == tmp_config_user_dir.join("sunpyrc").strpath


def test_get_and_create_download_dir(undo_download_dir_patch):
    # test default config
    path = get_and_create_download_dir()
    assert Path(path) == Path(USER) / 'sunpy' / 'data'
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
    assert Path(path) == Path(dirs.user_data_dir)
    # test updated config
    new_path = os.path.join(USER, 'sample_data_here_please')
    config.set('downloads', 'sample_dir', new_path)
    path = get_and_create_sample_dir()
    assert path == new_path
    # Set the config back
    os.rmdir(new_path)
    config.set('downloads', 'sample_dir', os.path.join(USER, 'sunpy', 'data', 'sample_data'))


def test_copy_default_config(tmpdir, undo_config_dir_patch, monkeypatch):
    config_filename = 'sunpyrc'

    # Try to set a manual one (already created)
    tmp_config_dir = tmpdir.mkdir("sunpy_test_configdir")
    monkeypatch.setenv("SUNPY_CONFIGDIR", tmp_config_dir.strpath)

    assert tmp_config_dir == _get_user_configdir()

    config_file = Path(sunpy.__file__).parent / 'data' / config_filename
    user_config_file = tmp_config_dir / config_filename

    # Create a new config file
    copy_default_config()
    with open(user_config_file) as ucf:
        with open(config_file) as cf:
            assert ucf.read() == cf.read()


def test_copy_default_config_without_overwrite(tmpdir, undo_config_dir_patch, monkeypatch):

    # Try to set a manual one (already created)
    tmp_config_dir = tmpdir.mkdir("sunpy_test_configdir")
    monkeypatch.setenv("SUNPY_CONFIGDIR", tmp_config_dir.strpath)

    assert tmp_config_dir == _get_user_configdir()

    # Create a new config file
    copy_default_config()

    # Without the `overwrite` parameter
    # the function should raise an error.
    with pytest.warns(SunpyUserWarning):
        copy_default_config()


def test_copy_default_config_with_overwrite(tmpdir, undo_config_dir_patch, monkeypatch):

    # Try to set a manual one (already created)
    tmp_config_dir = tmpdir.mkdir("sunpy_test_configdir")
    monkeypatch.setenv("SUNPY_CONFIGDIR", tmp_config_dir.strpath)

    assert tmp_config_dir == _get_user_configdir()

    # Create a new config file
    copy_default_config()

    # With the `overwrite` parameter
    # the function does not raise any error.
    with pytest.warns(SunpyUserWarning):
        copy_default_config(overwrite=True)


def test_copy_default_config_sitewarn(tmpdir, undo_config_dir_patch, monkeypatch):

    # Set a user config dir
    tmp_config_dir = tmpdir.mkdir("sunpy_test_configdir")
    monkeypatch.setenv("SUNPY_CONFIGDIR", tmp_config_dir.strpath)

    assert tmp_config_dir == _get_user_configdir()

    # Set a site config dir
    tmp_config_site_dir = tmpdir.mkdir("sunpy_test_configdir_site")
    monkeypatch.setattr(AppDirs, "site_config_dir", tmp_config_site_dir.strpath)

    # Create sunpyrc to site configdir
    tmp_file = tmp_config_site_dir.join("sunpyrc")
    # Have to write to the file otherwise its seen as a directory(?!)
    tmp_file.write("content")

    # With a sunpyrc file in the site config dir,
    # the function should raise a warning.
    with pytest.warns(SunpyUserWarning):
        copy_default_config()
