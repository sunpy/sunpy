"""
HyperMap
"""

from __future__ import absolute_import

import matplotlib.pyplot as plt

import sunpy.visualization as viz

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
        self.data = data.T
        self.system = coordinate_system
        self.header = header

    def plot(self, animate=False):
        """
        Plot a HyperMap

        Parameters
        ----------
        animate: int
            Axis of underlying ndarray to animate
        """
        #TODO: Make animate more cleverer
        naxis = len(self.system.frames)

        if naxis == 1:
            #Plot some thing
            pass
        elif naxis == 2:
            #imshow something
            pass
        elif naxis == 3:
            if isinstance(animate, bool) and animate:
                raise ValueError("What axis to animate fool?!")
            axlist = range(0, naxis)
            axlist.pop(animate)
            print axlist
            y_extent = self.system.frames[axlist[0]].get_extent()
            x_extent = self.system.frames[axlist[1]].get_extent()
            extent = x_extent + y_extent
            print extent
            #What axis are we animating over
            ax = plt.subplots(1)
            ani = viz.animate_array(self.data, animate, axes=ax, interval=200,
                                    cmap=plt.get_cmap('gray'), colorbar=False,
                                    norm='dynamic', extent=extent)

            return ani

        else:
            raise ValueError("Can't plot this hypermap, please slice it")
