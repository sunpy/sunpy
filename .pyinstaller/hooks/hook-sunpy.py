from PyInstaller.utils.hooks import collect_data_files, collect_entry_point, collect_submodules, copy_metadata

# Imports needed to run tests
datas, hiddenimports = collect_entry_point("pytest11")
hiddenimports += collect_submodules('numpy.distutils')
hiddenimports += collect_submodules('distutils')
hiddenimports += ['sunpy.data.data_manager.tests.mocks']
hiddenimports += collect_submodules('sunpy')
datas += collect_data_files("sunpy")
datas += collect_data_files("drms")
datas += copy_metadata("sunpy")
