"""
SJI IRIS
"""

from __future__ import absolute_import

__author__ = "Tomas Meszaros"
__email__ = "exo@tty.sk"

import sunpy.io
from sunpy.hypermap.coordinate_system import CoordinateFrame
from sunpy.hypermap.coordinate_system import CoordinateSystem
from sunpy.hypermap.coordinate_system import SpatialFrame

class Parser(object):
    """
    IRIS raw data parser.

    Parameters
    ----------
    iris_raw_data : string
        fits data file
    """

    def __init__(self, iris_raw_data):
        hdus = sunpy.io.read_file(iris_raw_data)
        self._header = hdus[0][1]

    _types_dictionary = {'HPLN-TAN': 'Spatial',
                         'HPLT-TAN': 'Spatial',
                         'Time': 'Time'}

    def _get_header_item_group(self, group):
        """
        Filter and return item group (e.g. 'CTYPE') from IRIS header items.
        """
        return filter(lambda x: not x[0].find(group) and
                                not x[0] == group, self._header.items())

    def _make_frame_list(self):
        """
        Make list of CoordinateFrames so that we can feed CoordinateSystem
        with it later.
        """
        frame_list = []
        num_axes = len(self._get_header_item_group('NAXIS'))
        ctypes = self._get_header_item_group('CTYPE')
        cunits = self._get_header_item_group('CUNIT')
        crpixs = self._get_header_item_group('CRPIX')
        # TODO: Where can we get axes names?

        for i in range(num_axes):
            frame_type = self._types_dictionary.get(ctypes[i][1])

            if frame_type == 'Spatial':
                frame_list.append(SpatialFrame(crpixs[i][1], None,
                                               [cunits[i][1], cunits[i][1]]))
            elif frame_type == 'Time':
                frame_list.append(CoordinateFrame('Time', 1, None,
                                                  [cunits[i][1]]))
        return frame_list

    def get_coordinate_system(self, name="CompositeSystem"):
        return CoordinateSystem(self._make_frame_list(), name)
