# Tomas Meszaros <exo@tty.sk>

from sunpy.hypermap.coordinate_system import CoordinateSystem
from sunpy.hypermap.coordinate_system import CoordinateFrame
from sunpy.hypermap.coordinate_system import SpatialFrame
from sunpy.hypermap.coordinate_system import SpectralFrame


# These tests only check creation of the *Frame objects, the attributes are
# checked further in the iris/test_iris.py tests.

def test_coordinate_frame():
    cframe_time = CoordinateFrame("time", 1.1, 2, 3, 4, ["time"], ["second"])
    assert isinstance(cframe_time, CoordinateFrame)

def test_spatial_frame():
    cframe_spatial = SpatialFrame(1.1, 2, 3, ["x", "y"], ["au", "au"])
    assert isinstance(cframe_spatial, SpatialFrame)

def test_spectral_frame():
    cframe_spectral = SpectralFrame(1.1, 2, 3, ["x", "y"], ["au", "au"])
    assert isinstance(cframe_spectral, SpectralFrame)

def test_coordinate_system():
    cframe_time = CoordinateFrame("time", 1.1, 2, 3, 4, ["time"], ["second"])
    cframe_spatial = SpatialFrame("Spatial", 1.1, 2, 3, ["x", "y"], ["au", "au"])
    cframe_spectral = SpectralFrame("Sectral", 1.1, 2, 3, ["whatever"], ["dummy"])
    frame_list = [cframe_time, cframe_spatial, cframe_spectral]
    system_name = "Random Test Coordinate System"
    csystem = CoordinateSystem(frame_list, system_name)
    assert isinstance(csystem, CoordinateSystem)
