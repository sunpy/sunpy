import sunpy

def test_manager_exists():
    assert hasattr(sunpy, "_data_manager")

def test_make_map_add():
    map0 = sunpy.make_map(sunpy.AIA_171_IMAGE)
    map1 = sunpy.make_map(sunpy.AIA_171_IMAGE)

    assert len(sunpy._data_manager.data) == 2
    assert isinstance(list(sunpy._data_manager.data)[0](), sunpy.Map)

def test_map_removed():
    map_ = sunpy.make_map(sunpy.AIA_171_IMAGE)
    assert len(sunpy._data_manager.data) == 1
    del map_
    assert len(sunpy._data_manager.data) == 0

