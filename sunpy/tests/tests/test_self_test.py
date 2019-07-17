import imp
import os.path
import warnings

import pytest

import sunpy
import sunpy.tests.runner
from sunpy.util.exceptions import SunpyDeprecationWarning

root_dir = os.path.dirname(os.path.abspath(sunpy.__file__))


def test_import_runner():
    """
    When running the tests with setup.py, the test runner class is imported by
    setup.py before coverage is watching.

    To ensure that the coverage for sunpy/tests/runner.py is correctly
    measured we force the interpreter to reload it here while coverage
    is watching.
    """
    imp.reload(sunpy.tests.runner)


def test_main_nonexisting_module():
    with pytest.raises(ValueError):
        sunpy.self_test(package='doesnotexist')


def test_main_stdlib_module():
    """
    This test makes sure that the module is really searched within the sunpy
    package.
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


def test_main_include_remote_data(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(package='map', online=True)
    assert args in ([os.path.join('sunpy', 'map'), '--remote-data=any', '-m', 'not figure'],
                    [os.path.join(root_dir, 'map'),'--remote-data=any',  '-m', 'not figure'])


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


def test_main_figure_dir(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(figure_only=True, figure_dir=".")
    assert args in (['sunpy', '--figure_dir', '.', '-m', 'figure'],
                    [root_dir, '--figure_dir', '.', '-m', 'figure'])


def test_main_coverage(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(coverage=True)
    for a in args:
        assert a in [root_dir, 'sunpy', '--cov', '--cov-config', '-m',
                     'not figure',
                     os.path.join(root_dir, 'tests', 'coveragerc'),
                     os.path.join('sunpy', 'tests', 'coveragerc')]


def test_main_coverage_report(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(coverage=True, cov_report=True)
    for a in args:
        assert a in [root_dir, 'sunpy', '--cov', '--cov-config', '-m',
                     'not figure',
                     os.path.join(root_dir, 'tests', 'coveragerc'),
                     os.path.join('sunpy', 'tests', 'coveragerc'),
                     '--cov-report']


def test_main_coverage_report_html(monkeypatch):
    monkeypatch.setattr(pytest, 'main', lambda args, **kwargs: args)
    args = sunpy.self_test(coverage=True, cov_report=True)
    for a in args:
        assert a in [root_dir, 'sunpy', '--cov', '--cov-config', '-m',
                     'not figure',
                     os.path.join(root_dir, 'tests', 'coveragerc'),
                     os.path.join('sunpy', 'tests', 'coveragerc'),
                     '--cov-report']


def test_warnings():
    # Ensure that our warning trickery dosen't stop pytest.warns working
    with pytest.warns(SunpyDeprecationWarning):
        warnings.warn("Hello", SunpyDeprecationWarning)
