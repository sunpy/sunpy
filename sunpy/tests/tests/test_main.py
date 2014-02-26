import os.path

import pytest

import sunpy.tests

root_dir = os.path.dirname(os.path.abspath(sunpy.__file__))

def test_main_nonexisting_module():
    with pytest.raises(ImportError):
        sunpy.tests.main('doesnotexist')


def test_main_stdlib_module():
    """This test makes sure that the module is really searched within the
    sunpy package.

    """
    with pytest.raises(ImportError):
        sunpy.tests.main('random')


def test_main_noargs(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda x: x)
    args = sunpy.tests.main()
    assert args == [root_dir]


def test_main_submodule(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda x: x)
    args = sunpy.tests.main('map')
    assert args == [os.path.join(root_dir, 'map', 'tests')]


def test_main_with_cover(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda x: x)
    args = sunpy.tests.main('map', cover=True)
    covpath = os.path.abspath(
        os.path.join(sunpy.tests.testdir, os.path.join(os.pardir, 'map')))
    assert args == ['--cov', covpath, os.path.join(root_dir, 'map', 'tests')]


def test_main_with_show_uncovered_lines(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda x: x)
    args = sunpy.tests.main('map', show_uncovered_lines=True)
    assert args == [
        '--cov-report', 'term-missing',
        os.path.join(root_dir, 'map', 'tests')]


def test_main_exclude_online(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda x: x)
    args = sunpy.tests.main('map', online=sunpy.tests.EXCLUDE_ONLINE)
    assert args == ['-k-online', os.path.join(root_dir, 'map', 'tests')]


def test_main_only_online(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda x: x)
    args = sunpy.tests.main('map', online=sunpy.tests.ONLY_ONLINE)
    assert args == ['-k', 'online', os.path.join(root_dir, 'map', 'tests')]


def test_main_invalid_online_parameter():
    with pytest.raises(ValueError) as excinfo:
        sunpy.tests.main(online='blabla')
    assert excinfo.exconly() == (
        'ValueError: `online` parameter must have one of the following '
        'values: sunpy.tests.INCLUDE_ONLINE, sunpy.tests.EXCLUDE_ONLINE, '
        'sunpy.tests.ONLY_ONLINE')
