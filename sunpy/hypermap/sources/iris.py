"""
SJI IRIS
"""

from __future__ import absolute_import

__author__ = "Tomas Meszaros"
__email__ = "exo@tty.sk"

import collections

from sunpy.hypermap.coordinate_system import CoordinateSystem
from sunpy.hypermap.coordinate_system import CoordinateFrame
from sunpy.hypermap.coordinate_system import SpatialFrame
from sunpy.hypermap.coordinate_system import SpectralFrame

class Parser(object):
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
            raise HeaderTypeError("Header have to be an instance of the \
                                   collections.Mapping!")
        self._header = header

    def _get_header_item_group(self, group):
        """
        Filter header and return list of items of a specific header
        group (e.g. 'CTYPE').
        Return empty list when unable to find @group in @_header.items().
        """
        return filter(lambda x: not x[0].find(group) and not x[0] == group,
                                self._header.items())

    def _make_frame_list(self):
        """
        Make list of CoordinateFrames so that we can feed CoordinateSystem
        with it later.
        """
        frame_list = []

        naxiss = self._get_header_item_group('NAXIS')
        ctypes = self._get_header_item_group('CTYPE')
        cunits = self._get_header_item_group('CUNIT')
        crpixs = self._get_header_item_group('CRPIX')
        cdelts = self._get_header_item_group('CDELT')
        # TODO: Where can we get axes names?

        # We want to stop when we are unable to find required keywords,
        # e.g. 'CTYPE', 'CUNIT', etc., in the header (header is probably
        # missing them).
        if [] in [naxiss, ctypes, cunits, crpixs]:
            raise MissingWCSKeywords("Missing WCS keyword in the header!")

        num_axes = len(self._get_header_item_group('NAXIS'))

        for i in range(num_axes):
            frame_type = self._types_dictionary.get(ctypes[i][1])

            if frame_type == 'Spatial':
                frame_list.append(SpatialFrame(reference_position=crpixs[i][1],
                                               pixel_size=cdelts[i][1],
                                               number_of_pixels=naxiss[i][1],
                                               axes_names=None,
                                               units=[cunits[i][1],
                                                      cunits[i][1]]))
            elif frame_type == 'Time':
                frame_list.append(CoordinateFrame(system='Time',
                                                  reference_position=crpixs[i][1],
                                                  pixel_size=cdelts[i][1],
                                                  number_of_pixels=naxiss[i][1],
                                                  num_axes=1,
                                                  axes_names=None,
                                                  units=[cunits[i][1]]))
            elif frame_type == 'Spectral':
                frame_list.append(SpectralFrame(reference_position=crpixs[i][1],
                                                pixel_size=cdelts[i][1],
                                                number_of_pixels=naxiss[i][1],
                                                axes_names=None,
                                                units=[cunits[i][1]]))
        return frame_list

    def get_coordinate_system(self, name="CompositeSystem"):
        return CoordinateSystem(self._make_frame_list(), name)


class HeaderTypeError(Exception):
    pass

class MissingWCSKeywords(Exception):
    pass
