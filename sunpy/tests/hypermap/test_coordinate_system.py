# Tomas Meszaros <exo@tty.sk>

from sunpy.hypermap.coordinate_system import CoordinateSystem
from sunpy.hypermap.coordinate_system import CoordinateFrame
from sunpy.hypermap.coordinate_system import SpatialFrame

def test_coordinate_frame():
    cframe_time = CoordinateFrame("time", 1, ["time"], ["second"])
    assert isinstance(cframe_time, CoordinateFrame)

def test_spatial_frame():
    cframe_spatial = SpatialFrame("Spatial", ["x", "y"], ["au", "au"])
    assert isinstance(cframe_spatial, SpatialFrame)

def test_coordinate_system():
    cframe_time = CoordinateFrame("time", 1, ["time"], ["second"])
    cframe_spatial = SpatialFrame("Spatial", ["x", "y"], ["au", "au"])
    frame_list = [cframe_time, cframe_spatial]
    system_name = "Random Test Coordinate System"
    csystem = CoordinateSystem(frame_list, system_name)
    assert isinstance(csystem, CoordinateSystem)
