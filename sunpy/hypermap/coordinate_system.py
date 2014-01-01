""" Coordinate system for HyperMap. """
from __future__ import absolute_import

__author__ = "Tomas Meszaros"
__email__ = "exo@tty.sk"

import collections


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

    def __init__(self, system, reference_position, reference_coordinate,
                 pixel_size, number_of_pixels, num_axes, axes_names=None,
                 units=None):
        """ Initialize a frame"""
        self.system = system
        self.reference_position = reference_position
        self.reference_coordinate = reference_coordinate
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

    def get_extent(self):
        """
        Return a list of the min and max values for the axis
        """
        amin = (self.reference_position - self.number_of_pixels / 2. *
                self.pixel_size)
        amax = (self.reference_position + self.number_of_pixels / 2. *
                self.pixel_size)
        return [amin, amax]


class SpatialFrame(CoordinateFrame):
    """
    SpatialFrame

    Parameters
    -----------------
    reference_position: list, BUT possible astropy.Time or coordinate object
                        in the future
    """

    def __init__(self, reference_position, reference_coordinate, pixel_size,
                 number_of_pixels, axes_names=["",""], units=["",""]):
        super(SpatialFrame, self).__init__(
          system='Spatial',
          reference_position=reference_position,
          reference_coordinate=reference_coordinate,
          pixel_size=pixel_size,
          number_of_pixels=number_of_pixels,
          num_axes=2,
          axes_names=axes_names,
          units=units
        )


class SpectralFrame(CoordinateFrame):
    """
    SpectralFrame

    Parameters
    -----------------
    reference_position: list, BUT possible astropy.Time or coordinate object
                        in the future
    """

    def __init__(self, reference_position, reference_coordinate, pixel_size,
                 number_of_pixels, axes_names=["",""], units=["",""]):
        super(SpectralFrame, self).__init__(
          system='Spectral',
          reference_position=reference_position,
          reference_coordinate=reference_coordinate,
          pixel_size=pixel_size,
          number_of_pixels=number_of_pixels,
          num_axes=1,
          axes_names=axes_names,
          units=units
        )


class WCSParser(object):
    """
    IRIS raw data parser.

    Parameters
    ----------
    header : collections.Mapping instance
        e.g.: sunpy.io.header.FileHeader
    """

    _types_dictionary = {'HPLN-TAN': 'Spatial',
                         'HPLT-TAN': 'Spatial',
                         'Time': 'Time',
                         'WAVE': 'Spectral'}

    def __init__(self, header):
        if not isinstance(header, collections.Mapping):
            raise HeaderTypeError("Header have to be an instance of the "
                                  "collections.Mapping!")
        self._header = header

    def get_header_item_group(self, group):
        """
        Filter header and return list of items of a specific header
        group (e.g. 'CTYPE').
        Return empty list when unable to find @group in @_header.items().
        """
        return [i for i in self._header.items() if not i[0].find(group) and
                                                   not i[0] == group]

    def _make_frame_list(self, reverse=False):
        """
        Make list of CoordinateFrames so that we can feed CoordinateSystem
        with it later.
        """
        frames = []

        naxiss = self.get_header_item_group('NAXIS')
        ctypes = self.get_header_item_group('CTYPE')
        cunits = self.get_header_item_group('CUNIT')
        crpixs = self.get_header_item_group('CRPIX')
        cdelts = self.get_header_item_group('CDELT')
        crvals = self.get_header_item_group('CRVAL')

        # TODO: Where can we get axes names?

        # We want to stop when we are unable to find required keywords,
        # e.g. 'CTYPE', 'CUNIT', etc., in the header (header is probably
        # missing them).
        if [] in [naxiss, ctypes, cunits, crpixs]:
            raise MissingWCSKeywords("Missing WCS keyword in the header!")

        num_axes = len(self.get_header_item_group('NAXIS'))

        for i in range(num_axes):
            frame_type = self._types_dictionary.get(ctypes[i][1])

            if frame_type == 'Spatial':
                frames.append(SpatialFrame(reference_position=crpixs[i][1],
                                           reference_coordinate=crvals[i][1],
                                           pixel_size=cdelts[i][1],
                                           number_of_pixels=naxiss[i][1],
                                           axes_names=None,
                                           units=[cunits[i][1],
                                                  cunits[i][1]]))
            elif frame_type == 'Time':
                frames.append(CoordinateFrame(system='Time',
                                              reference_position=crpixs[i][1],
                                              reference_coordinate=crvals[i][1],
                                              pixel_size=cdelts[i][1],
                                              number_of_pixels=naxiss[i][1],
                                              num_axes=1,
                                              axes_names=None,
                                              units=[cunits[i][1]]))
            elif frame_type == 'Spectral':
                frames.append(SpectralFrame(reference_position=crpixs[i][1],
                                            reference_coordinate=crvals[i][1],
                                            pixel_size=cdelts[i][1],
                                            number_of_pixels=naxiss[i][1],
                                            axes_names=None,
                                            units=[cunits[i][1]]))
        if reverse:
            return frames[::-1]
        else:
            return frames

    def get_coordinate_system(self, name="CompositeSystem", reverse=True):
        """
        Create CoordinateSystem out of frame list.

        Parameters
        ----------
        name: str
            Name for the created system

        reverse: bool
            Reverse the order of the frames, WCS headers are by default transposed
            to their numpy array axis ordering.
        """
        return CoordinateSystem(self._make_frame_list(reverse=reverse), name)


class HeaderTypeError(Exception):
    """ Exception to raise when header is not an instance of collections.Mapping
    """
    pass


class MissingWCSKeywords(Exception):
    """ Exception to raise when is header missing some WCS keyword.
    """
    pass
