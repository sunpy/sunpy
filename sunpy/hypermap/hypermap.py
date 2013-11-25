"""
HyperMap
"""

from __future__ import absolute_import

__author__ = "Tomas Meszaros"
__email__ = "exo@tty.sk"

class HyperMap(object):
    """
    HyperMap

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
