"""
HyperMap
"""
from __future__ import absolute_import

__author__ = "Tomas Meszaros"
__email__ = "exo@tty.sk"

from sunpy.hypermap.coordinate_system import CoordinateSystem
from sunpy.hypermap.coordinate_system import CoordinateFrame


class Hypermap(object):
    """
    Hypermap

    Parameters
    ----------
    data : numpy.ndarray
        data
    coordinate_system : CoordinateSystem
        CoordinateSystem object
    header : list
        header
    """
    def __init__(self, data, coordinate_system, header):
        self.data = data
        self.system = coordinate_system
        self.header = header

