from PyInstaller.utils.hooks import collect_data_files, collect_entry_point, collect_submodules, copy_metadata

datas, hiddenimports = collect_entry_point("pytest11")

datas += collect_data_files("sunpy")
datas += collect_data_files("dask")
datas += collect_data_files("drms")
datas += copy_metadata("sunpy")

hiddenimports += collect_submodules('sunpy')
hiddenimports += collect_submodules('numpy.distutils')
hiddenimports += collect_submodules('distutils')
# hiddenimports += ['skimage.filters.rank.core_cy_3d']
hiddenimports += ['sunpy.data.data_manager.tests.mocks']
