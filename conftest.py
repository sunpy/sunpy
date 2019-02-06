import importlib

from pkg_resources import parse_version

if importlib.util.find_spec('asdf') is not None:
    from asdf import __version__ as asdf_version
    if parse_version(asdf_version) >= parse_version('2.3.0'):
        pytest_plugins = ['asdf.tests.schema_tester']
