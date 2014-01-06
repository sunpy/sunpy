"""
HyperMap
"""

from __future__ import absolute_import

import matplotlib.pyplot as plt

import sunpy.visualization as viz
from sunpy.visualization.imageanimator import ImageAnimator

from . coordinate_system import WCSParser, SpatialFrame

__author__ = "Tomas Meszaros"
__email__ = "exo@tty.sk"


class HyperMap(object):
    """
    HyperMap

    Parameters
    ----------
    data : numpy.ndarray
        data
    header : list
        header
    coordinate_system : CoordinateSystem
        CoordinateSystem object
    """

    def __init__(self, data, header, coordinate_system=False):
        self.data = data
        self.header = header

        if not coordinate_system:
            self.system = WCSParser(header).get_coordinate_system()
        else:
            self.system = coordinate_system

    def plot(self, animate=False, **kwargs):
        """
        Plot a HyperMap

        This function calls many different underlying plotting routines
        depending on the data in the hypermap.

        Parameters
        ----------
        animate: int
            Axis of underlying ndarray to animate
        """
        naxis = len(self.system.frames)
        num_spatial = sum([isinstance(frame, SpatialFrame) for frame in self.system.frames])

        if num_spatial == 2 and naxis > 2 and not animate:
            ImageAnimator(self.data, axis_range=[s.get_extent() for s in self.system.frames], **kwargs)
        else:
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
                y_extent = self.system.frames[axlist[0]].get_extent()
                x_extent = self.system.frames[axlist[1]].get_extent()
                extent = x_extent + y_extent
                #setup kwargs
                anim_kwargs = {'interval':200,
                               'colorbar':True,
                               'norm':'dynamic',
                               'extent':extent}
                anim_kwargs.update(kwargs)
                #What axis are we animating over
                ax = plt.subplot()
                ani = viz.animate_array(self.data, animate, axes=ax, **anim_kwargs)

                return ani

            else:
                raise ValueError("Can't plot this hypermap, please slice it")
