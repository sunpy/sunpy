import os.path

import pytest

import sunpy

root_dir = os.path.dirname(os.path.abspath(sunpy.__file__))


def test_main_nonexisting_module():
    with pytest.raises(ValueError):
        sunpy.self_test(package='doesnotexist')


def test_main_stdlib_module():
    """This test makes sure that the module is really searched within the
    sunpy package.

    """
    with pytest.raises(ValueError):
        sunpy.self_test(package='random')


def test_main_noargs(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test()
    assert args in (['sunpy', '-m', 'not figure'],
                    [root_dir, '-m', 'not figure'])


def test_main_submodule_map(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(package='map')
    assert args in ([os.path.join('sunpy', 'map'), '-m', 'not figure'],
                    [os.path.join(root_dir, 'map'), '-m', 'not figure'])


def test_main_submodule_jsoc(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(package='net.jsoc')
    assert args in ([os.path.join('sunpy', 'net', 'jsoc'), '-m', 'not figure'],
                    [os.path.join(root_dir, 'net', 'jsoc'), '-m', 'not figure'])


def test_main_exclude_remote_data(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(package='map', online=False)
    assert args in ([os.path.join('sunpy', 'map'), '-m', 'not figure'],
                    [os.path.join(root_dir, 'map'), '-m', 'not figure'])


def test_main_only_remote_data(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(package='map', online_only=True)
    assert args in ([os.path.join('sunpy', 'map'), '-k remote_data', '--remote-data=any', '-m', 'not figure'],
                    [os.path.join(root_dir, 'map'), '-k remote_data', '--remote-data=any', '-m', 'not figure'])


def test_main_figures(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(figure=True)
    assert args in (['sunpy'],
                    [root_dir])


def test_main_figure_only(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(figure_only=True)
    assert args in (['sunpy', '-m', 'figure'],
                    [root_dir, '-m', 'figure'])
