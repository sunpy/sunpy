from packaging.requirements import Requirement

from sunpy.util.sysinfo import (
    find_dependencies,
    format_requirement_string,
    missing_dependencies_by_extra,
    resolve_requirement_versions,
    system_info,
)

EXTRA_DEPS = [
    'asdf-astropy',
    'asdf',
    'astropy',
    'beautifulsoup4',
    'cdflib',
    'drms',
    'h5netcdf',
    'h5py',
    'lxml',
    'matplotlib',
    'mpl-animators',
    'numpy',
    'packaging',
    'pandas',
    'parfive',
    'pyerfa',
    'python-dateutil',
    'reproject',
    'scipy',
    'zeep',
]

EXTRA_ALL_GROUPS = [
    'all',
    'asdf',
    'dev',
    'docs-gallery',
    'docs',
    'image',
    'jpeg2000',
    'map',
    'net',
    'opencv',
    'required',
    'spice',
    'tests',
    'timeseries',
    'visualization',
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

    # There should not be any self-referential dependency on sunpy
    assert "sunpy" not in installed


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
    lines = captured.out.splitlines()
    assert "sunpy Installation Information" in lines

    # sunpy should not be listed as an optional dependency
    index = lines.index("Optional Dependencies")
    for line in lines[index + 2:]:
        assert not line.startswith("sunpy:")
