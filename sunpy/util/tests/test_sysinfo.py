from packaging.requirements import Requirement

from sunpy.util.sysinfo import (
    find_dependencies,
    format_requirement_string,
    missing_dependencies_by_extra,
    resolve_requirement_versions,
    system_info,
)

EXTRA_DEPS = [
    'asdf',
    'asdf-astropy',
    'astropy',
    'numpy',
    'parfive',
    'packaging',
    'dask',
    'scikit-image',
    'scipy',
    'lxml',
    'matplotlib',
    'mpl-animators',
    'reproject',
    'beautifulsoup4',
    'drms',
    'python-dateutil',
    'zeep',
    'cdflib',
    'h5netcdf',
    'h5py',
    'pandas'
]

EXTRA_ALL_GROUPS = [
    'all',
    'asdf',
    'required',
    'dask',
    'dev',
    'docs',
    'docs-gallery',
    'image',
    'jpeg2000',
    'map',
    'net',
    'spice',
    'tests',
    'timeseries',
    'visualization'
]

def test_find_dependencies():
    """
    This is ran in several test environments with varying dependencies installed.
    So it will be common to find not docs installed, so there will be "missing" dependencies.
    But this is not a problem.
    """
    _, installed = find_dependencies(package="sunpy", extras=["required", *EXTRA_ALL_GROUPS])
    for package in EXTRA_DEPS:
        assert package in installed


def test_missing_dependencies_by_extra():
    missing = missing_dependencies_by_extra()
    for group in EXTRA_ALL_GROUPS:
        assert group in missing


def test_resolve_requirement_versions():
    package1 = Requirement('test-package[ext1]>=1.1.1; extra == "group1"')
    package2 = Requirement('test-package[ext2]<=2.0.0; extra == "group2"')
    assert str(resolve_requirement_versions([package1, package2])) == str(Requirement(
        'test-package[ext1,ext2]<=2.0.0,>=1.1.1; extra == "group1" or extra == "group2"'))

    package3 = Requirement('test-package==1.1.0; extra == "group3"')
    package4 = Requirement('test-package==1.1.0; extra == "group4"')
    assert str(resolve_requirement_versions([package3, package4])) == str(
        Requirement('test-package==1.1.0; extra == "group3" or extra == "group4"'))

    package5 = Requirement('test-package; extra == "group5"')
    package6 = Requirement('test-package[ext3]@https://foo.com')
    assert str(resolve_requirement_versions([package5, package6])) == str(
        Requirement('test-package[ext3]@ https://foo.com ; extra == "group5"'))


def test_format_requirement_string():
    package1 = Requirement('test-package[ext1]>=1.1.1; extra == "group1"')
    assert format_requirement_string(package1) == 'Missing test-package[ext1]>=1.1.1; extra == "group1"'

    package2 = Requirement('test-package>=1.1.1; extra == "group1" or extra == "group2" or extra == "group3"')
    assert format_requirement_string(
        package2) == 'Missing test-package>=1.1.1; extra == "group1" or "group2" or "group3"'

    package3 = Requirement('test-package>=1.1.1')
    assert format_requirement_string(package3) == 'Missing test-package>=1.1.1'


def test_system_info(capsys):
    system_info()
    captured = capsys.readouterr()
    assert "\nsunpy Installation Information\n" in captured.out
