from sunpy.util.sysinfo import find_dependencies, missing_dependencies_by_extra, system_info


def test_find_dependencies():
    missing, installed = find_dependencies()
    assert missing == {}
    assert sorted(list(installed.keys())) == sorted(["astropy", "numpy", "packaging", "parfive"])

    missing, installed = find_dependencies(package="sunpy", extras=["required", "all"])
    assert missing == {}
    assert sorted(list(installed.keys())) == sorted(['astropy', 'numpy', 'parfive', 'packaging', 'asdf',
                                                     'dask', 'sqlalchemy', 'scikit-image', 'scipy', 'glymur',
                                                     'matplotlib', 'mpl-animators', 'reproject',
                                                     'beautifulsoup4', 'drms', 'python-dateutil', 'tqdm',
                                                     'zeep', 'cdflib', 'h5netcdf', 'h5py', 'pandas'])


def test_missing_dependencies_by_extra():
    missing = missing_dependencies_by_extra()
    assert sorted(list(missing.keys())) == sorted(['all', 'asdf', 'required', 'dask', 'database', 'dev', 'docs',
                                                   'image', 'jpeg2000', 'map', 'net', 'tests', 'timeseries',
                                                   'visualization'])
    missing = missing_dependencies_by_extra(exclude_extras=["all"])
    assert sorted(list(missing.keys())) == sorted(['asdf', 'required', 'dask', 'database', 'dev', 'docs',
                                                   'image', 'jpeg2000', 'map', 'net', 'tests', 'timeseries',
                                                   'visualization'])


def test_system_info(capsys):
    system_info()
    captured = capsys.readouterr()
    assert "\nsunpy Installation Information\n" in captured.out
