"""
Coordinate system for HyperMap.
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
    reference_position : float
    pixel_size : float
    number_of_pixels : int
    num_axes : int
    axes_names : list of strings
    units : list of units
    """

    def __init__(self, system, reference_position, pixel_size, number_of_pixels,
                 num_axes, axes_names=None, units=None):
        """ Initialize a frame"""
        self.system = system
        self.reference_position = reference_position
        self.pixel_size = pixel_size
        self.number_of_pixels = number_of_pixels
        self.num_axes = num_axes
        self.axes_names = axes_names
        self.units = units

    def transform_to(self, other):
        """
        Transform from the current reference system to other if
        the system attribute of the two matches
        """


class SpatialFrame(CoordinateFrame):
    """
    SpatialFrame

    Parameters
    -----------------
    reference_position: list, BUT possible astropy.Time or coordinate object
                        in the future
    """

    def __init__(self, reference_position, pixel_size, number_of_pixels,
                 axes_names=["",""], units=["",""]):
        super(SpatialFrame, self).__init__(system='Spatial',
                                           reference_position=reference_position,
                                           pixel_size=pixel_size,
                                           number_of_pixels=number_of_pixels,
                                           num_axes=2,
                                           axes_names=axes_names,
                                           units=units)


class SpectralFrame(CoordinateFrame):
    """
    SpectralFrame

    Parameters
    -----------------
    reference_position: list, BUT possible astropy.Time or coordinate object
                        in the future
    """

    def __init__(self, reference_position, pixel_size, number_of_pixels,
                 axes_names=["",""], units=["",""]):
        super(SpectralFrame, self).__init__(system='Spectral',
                                           reference_position=reference_position,
                                           pixel_size=pixel_size,
                                           number_of_pixels=number_of_pixels,
                                           num_axes=1,
                                           axes_names=axes_names,
                                           units=units)
