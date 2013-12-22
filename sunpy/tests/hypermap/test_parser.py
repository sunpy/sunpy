# Tomas Meszaros <exo@tty.sk>

""" Tests for hypermap/coordinate_system.py WCSParser
"""

# NOTE: All test are "disabled" until the test data will be available.

import sunpy.io  # read_file blows away RAM when using big Raster fits
from sunpy.hypermap.coordinate_system import CoordinateSystem, CoordinateFrame
from sunpy.hypermap.coordinate_system import WCSParser


# This is temporary, Stuart Mumford will provide nice & small test samples :-)
_IRIS_SAMPLE_DATA = "/home/exo/iris_sample_data.fits"
_IRIS_RASTER_SAMPLE_DATA = "/home/exo/iris_raster_sample_data.fits"


def test_item_group():
    """ Test for correct header item group extraction.
    """
    hdus = sunpy.io.read_file(_IRIS_SAMPLE_DATA)
    header = hdus[0][1]
    parsed = WCSParser(header)

    assert parsed.get_header_item_group('DUMMY_FOO_JANE_DOE_#!@') == []
    assert parsed.get_header_item_group('NAXIS') != []
    assert parsed.get_header_item_group('CTYPE') != []
    assert parsed.get_header_item_group('CUNIT') != []
    assert parsed.get_header_item_group('CRPIX') != []
    assert parsed.get_header_item_group('CDELT') != []


def test_raster_item_group():
    """ Test for correct raster header item group extraction.
    """
    headers = sunpy.io.read_file_header(_IRIS_RASTER_SAMPLE_DATA)
    parsed_0 = WCSParser(headers[0])

    assert parsed_0.get_header_item_group('DUMMY_FOO_JANE_DOE_#!@') == []
    assert parsed_0.get_header_item_group('NAXIS') == []
    assert parsed_0.get_header_item_group('CTYPE') == []
    assert parsed_0.get_header_item_group('CUNIT') == []
    assert parsed_0.get_header_item_group('CRPIX') == []
    assert parsed_0.get_header_item_group('CDELT') == []

    parsed_1 = WCSParser(headers[1])

    assert parsed_1.get_header_item_group('DUMMY_FOO_JANE_DOE_#!@') == []
    assert parsed_1.get_header_item_group('NAXIS') != []
    assert parsed_1.get_header_item_group('CTYPE') != []
    assert parsed_1.get_header_item_group('CUNIT') != []
    assert parsed_1.get_header_item_group('CRPIX') != []
    assert parsed_1.get_header_item_group('CDELT') != []


def test_make_coord_system():
    """ Test CoordinateSystem creation and functionality.
    """
    hdus = sunpy.io.read_file(_IRIS_SAMPLE_DATA)
    header = hdus[0][1]
    parsed = WCSParser(header)
    system = parsed.get_coordinate_system("Coordinate Test System")

    assert isinstance(system, CoordinateSystem)
    for i in system.frames:
        assert isinstance(i, CoordinateFrame)
        assert type(i.system) == str and i.system != ""
        assert type(i.pixel_size) == float and i.pixel_size >= 0
        assert type(i.number_of_pixels) == int and i.number_of_pixels > 0
        assert type(i.num_axes) == int and i.num_axes > 0
        assert type(i.axes_names) == list or i.axes_names is None
        assert type(i.units) == list or i.units is None


def test_make_raster_coord_system():
    """ Test raster CoordinateSystem creation and functionality.
    """
    headers = sunpy.io.read_file_header(_IRIS_RASTER_SAMPLE_DATA)
    parsed = WCSParser(headers[1])
    system = parsed.get_coordinate_system("Coordinate Test System")

    assert isinstance(system, CoordinateSystem)
    for i in system.frames:
        assert isinstance(i, CoordinateFrame)
        assert type(i.system) == str and i.system != ""
        assert type(i.pixel_size) == float and i.pixel_size > 0
        assert type(i.number_of_pixels) == int and i.number_of_pixels > 0
        assert type(i.num_axes) == int and i.num_axes > 0
        assert type(i.axes_names) == list or i.axes_names is None
        assert type(i.units) == list or i.units is None


# ================== #
# Just for debugging #
# ================== #

def _show_coordinate_system():
    """ Print created CoordinateSystem.
    """
    hdus = sunpy.io.read_file(_IRIS_SAMPLE_DATA)
    header = hdus[0][1]
    parsed = WCSParser(header)
    system = parsed.get_coordinate_system("Coordinate Test System")

    print "=" * len(system.name)
    print system.name
    print "=" * len(system.name)

    for frame in system.frames:
        print "            system: %s" % frame.system
        print "reference_position: %f" % frame.reference_position
        print "        pixel_size: %f" % frame.pixel_size
        print "  number_of_pixels: %d" % frame.number_of_pixels
        print "          num_axes: %d" % frame.num_axes
        print "        axes_names: %s" % frame.axes_names
        print "             units: %s" % frame.units
        print "=" * len(system.name)


def _show_raster_coordinate_system():
    """ Print created raster CoordinateSystem.
    """
    headers = sunpy.io.read_file_header(_IRIS_RASTER_SAMPLE_DATA)

    for i in range(1, 5):
        parsed = WCSParser(headers[i])
        system = parsed.get_coordinate_system("Coordinate System " + str(i))

        print "=" * len(system.name)
        print system.name
        print "=" * len(system.name)

        for frame in system.frames:
            print "            system: %s" % frame.system
            print "reference_position: %f" % frame.reference_position
            print "        pixel_size: %f" % frame.pixel_size
            print "  number_of_pixels: %d" % frame.number_of_pixels
            print "          num_axes: %d" % frame.num_axes
            print "        axes_names: %s" % frame.axes_names
            print "             units: %s" % frame.units
            print "=" * len(system.name)
