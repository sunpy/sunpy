"""
Coordinate system for HyperMap
"""
from __future__ import absolute_import

__author__ = "Tomas Meszaros"
__email__ = "exo@tty.sk"


class CoordinateSystem(object):
    """
    A coordinate system has one or more frames.

    Parameters
    ----------
    name : string
        a user defined name
    frames : list
        list of CoordinateFrames [Time, Sky, Spectral]
    """
    def __init__(self, frames, name="CompositeSystem"):
        self.frames = frames
        self.name = name


class CoordinateFrame(object):
    """
    Base class for CoordinateFrames

    Parameters
    ----------
    system : string
        type of the frame
    num_axes : int
    axes_names : list of strings
    units : list of units
    """

    def __init__(self, system, num_axes, axes_names=None, units=None):
        """ Initialize a frame"""
        # TODO: Write type checking for @system, @num_axes, @axes_names, @units!
        self.system = system
        self.num_axes = num_axes
        self.axes_names = axes_names
        self.units = units

    def transform_to(self, other):
        """
        Transform from the current reference system to other if
        the system attribute of the two matches
        """


class SpatialFrame(CoordinateFrame):

    def __init__(self, reference_position, axes_names=["", ""], units=["",""]):
        """
        Parameters
        -----------------
        reference_position: list, BUT possible astropy.Time or coordinate
                                  object in the future
        """

        super(SpatialFrame, self).__init__('spatial', num_axes=2, axes_names=axes_names, units=units)
        self._reference_position = reference_position


# TESTS #
def run_test():
    # Create frames.
    cframe_time = CoordinateFrame("time", 1, ["time"], ["second"])
    cframe_spatial = SpatialFrame([0, 0], ["x", "y"], ["au", "au"])

    # Create coordinate system out of frames
    frame_list = [cframe_time, cframe_spatial]
    system_name = "Random Test Coordinate System"
    csystem = CoordinateSystem(frame_list, system_name)

    # Print stuff.
    test_print_coordinate_system(csystem)

def test_print_coordinate_system(csystem):
    print " Coordinate system: %s\n" \
          "  Number of frames: %s\n" \
          " -----------------\n" % (csystem.name, len(csystem.frames))

    for frame in csystem.frames:
        test_print_frame(frame)

def test_print_frame(cframe):
    format_touple = (cframe.system, cframe.num_axes, cframe.axes_names,
                     cframe.units)
    print "            system: %s\n" \
          "          num_axes: %d\n" \
          "        axes_names: %s\n" \
          "             units: %s" % format_touple
    if isinstance(cframe, SpatialFrame):
        print "reference_position: %s" % cframe._reference_position
    print "\n"
