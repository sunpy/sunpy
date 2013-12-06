# Tomas Meszaros <exo@tty.sk>

# NOTE: All test are "disabled" until the test data will be available.

import sunpy.io  #read_file blows away RAM when using big Raster fits
from sunpy.hypermap.sources.iris import Parser
from sunpy.hypermap.coordinate_system import CoordinateSystem
from sunpy.hypermap.coordinate_system import CoordinateFrame

# This is temporary, Stuart Mumford will provide nice & small test samples :-)
_IRIS_SAMPLE_DATA = "/home/exo/iris_sample_data.fits"
_IRIS_RASTER_SAMPLE_DATA = "/home/exo/iris_raster_sample_data.fits"

def _test_get_header_item_group():
    hdus = sunpy.io.read_file(_IRIS_SAMPLE_DATA)
    header = hdus[0][1]
    p = Parser(header)

    assert p._get_header_item_group('DUMMY_FOO_JANE_DOE_#!@') == []
    assert p._get_header_item_group('NAXIS') != []
    assert p._get_header_item_group('CTYPE') != []
    assert p._get_header_item_group('CUNIT') != []
    assert p._get_header_item_group('CRPIX') != []

def _test_get_header_item_group_raster():
    headers = sunpy.io.read_file_header(_IRIS_RASTER_SAMPLE_DATA)
    p0 = Parser(headers[0])

    assert p0._get_header_item_group('DUMMY_FOO_JANE_DOE_#!@') == []
    assert p0._get_header_item_group('NAXIS') == []
    assert p0._get_header_item_group('CTYPE') == []
    assert p0._get_header_item_group('CUNIT') == []
    assert p0._get_header_item_group('CRPIX') == []

    p1 = Parser(headers[1])

    assert p1._get_header_item_group('DUMMY_FOO_JANE_DOE_#!@') == []
    assert p1._get_header_item_group('NAXIS') != []
    assert p1._get_header_item_group('CTYPE') != []
    assert p1._get_header_item_group('CUNIT') != []
    assert p1._get_header_item_group('CRPIX') != []

def _test_make_coord_system():
    hdus = sunpy.io.read_file(_IRIS_SAMPLE_DATA)
    header = hdus[0][1]
    p = Parser(header)
    s = p.get_coordinate_system("Coordinate Test System")

    assert isinstance(s, CoordinateSystem)
    for i in s.frames:
        assert isinstance(i, CoordinateFrame)
        assert type(i.system) == str and i.system != ""
        assert type(i.num_axes) == int and i.num_axes > 0
        assert type(i.axes_names) == list or i.axes_names is None
        assert type(i.units) == list or i.units is None

def _test_make_raster_coord_system():
    headers = sunpy.io.read_file_header(_IRIS_RASTER_SAMPLE_DATA)
    p1 = Parser(headers[1])
    s1 = p1.get_coordinate_system("Coordinate Test System")

    assert isinstance(s1, CoordinateSystem)
    for i in s1.frames:
        assert isinstance(i, CoordinateFrame)
        assert type(i.system) == str and i.system != ""
        assert type(i.num_axes) == int and i.num_axes > 0
        assert type(i.axes_names) == list or i.axes_names is None
        assert type(i.units) == list or i.units is None
