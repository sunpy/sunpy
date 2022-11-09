
import pytest

import sunpy
from sunpy.tests.self_test import _self_test_args
from sunpy.util.exceptions import SunpyDeprecationWarning, warn_deprecated


def test_main_nonexisting_module():
    with pytest.raises(ModuleNotFoundError):
        sunpy.self_test(package='doesnotexist')


def test_main_stdlib_module():
    """
    This test makes sure that the module is really searched within the sunpy
    package.
    """
    with pytest.raises(ModuleNotFoundError):
        sunpy.self_test(package='random')


def test_main_noargs(monkeypatch):
    test_args = _self_test_args()
    assert test_args == ['-W', 'ignore', '--pyargs', 'sunpy']


def test_main_submodule_map(monkeypatch):
    args = _self_test_args(package='map')
    assert args == ['-W', 'ignore', '--pyargs', 'sunpy.map']


def test_main_submodule_jsoc(monkeypatch):
    args = _self_test_args(package='net.jsoc')
    assert args == ['-W', 'ignore', '--pyargs', 'sunpy.net.jsoc']


def test_main_exclude_remote_data(monkeypatch):
    args = _self_test_args(package='map', online=False)
    assert args == ['-W', 'ignore', '--pyargs', 'sunpy.map']


def test_main_include_remote_data(monkeypatch):
    args = _self_test_args(package='map', online=True)
    assert args == ['-W', 'ignore', '--remote-data=any', '--pyargs', 'sunpy.map']


def test_main_only_remote_data(monkeypatch):
    args = _self_test_args(package='map', online_only=True)
    assert args == ['-W', 'ignore', '--remote-data=any -m remote_data', '--pyargs', 'sunpy.map']


def test_main_figure_only(monkeypatch):
    args = _self_test_args(figure_only=True)
    assert args == ['-W', 'ignore', '--pyargs', 'sunpy', '-m', 'mpl_image_compare']


def test_warnings():
    # Ensure that our warning trickery doesn't stop pytest.warns working
    with pytest.warns(SunpyDeprecationWarning):
        warn_deprecated("Hello")
