# Tomas Meszaros <exo@tty.sk>

""" Tests for hypermap/coorinate_system.py
"""

from sunpy.hypermap.coordinate_system import CoordinateSystem, CoordinateFrame
from sunpy.hypermap.coordinate_system import SpatialFrame, SpectralFrame

# These tests only check creation of the *Frame objects, the attributes are
# checked further in the iris/test_iris.py tests.


def test_coordinate_frame():
    """ Test CoordinateFrame
    """
    cframe_time = CoordinateFrame("time", 1.1, 2, 3, 4, ["time"], ["second"])
    assert isinstance(cframe_time, CoordinateFrame)


def test_spatial_frame():
    """ Test SpatialFrame
    """
    cframe_spatial = SpatialFrame(1.1, 2, 3, ["x", "y"], ["au", "au"])
    assert isinstance(cframe_spatial, SpatialFrame)


def test_spectral_frame():
    """ Test SpectralFrame
    """
    cframe_spectral = SpectralFrame(1.1, 2, 3, ["x", "y"], ["au", "au"])
    assert isinstance(cframe_spectral, SpectralFrame)


def test_coordinate_system():
    """ Test CoordinateSystem
    """
    time = CoordinateFrame("time", 1.1, 2, 3, 4, ["time"], ["second"])
    spatial = SpatialFrame("Spatial", 1.1, 2, 3, ["x", "y"], ["au", "au"])
    spectral = SpectralFrame("Sectral", 1.1, 2, 3, ["whatever"], ["dummy"])
    frame_list = [time, spatial, spectral]
    system_name = "Random Test Coordinate System"
    csystem = CoordinateSystem(frame_list, system_name)
    assert isinstance(csystem, CoordinateSystem)
