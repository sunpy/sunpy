"""A Composite Map class

Author: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
from __future__ import absolute_import

import sunpy
import matplotlib.pyplot as plt
from sunpy.map.basemap import BaseMap
from sunpy.map.sources.rhessi import RHESSIMap

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

class CompositeMap:
    """
    CompositeMap(map1 [,map2,..])
    
    Parameters
    ----------
    args : *{sunpy.map, string}
        One or more map of filepaths
    
    Examples
    --------
    >>> import sunpy
    >>> sunpy.CompositeMap(sunpy.AIA_171_IMAGE, sunpy.RHESSI_IMAGE).show()
        
    >>> comp_map = sunpy.CompositeMap(sunpy.AIA_171_IMAGE, sunpy.EIT_195_IMAGE)    
    >>> comp_map.add_map(sunpy.RHESSI_IMAGE)
    >>> comp_map.show()
    
    """    
    def __init__(self, *args):
        self._maps = []
        
        # Default alpha and zorder values
        alphas = [1] * len(args)
        zorders = range(0, 10 * len(args), 10) 
        
        # Parse input Maps/filepaths        
        for i, item in enumerate(args):
            # Parse map
            if isinstance(item, BaseMap):
                m = item
            else:
                m = BaseMap.read(item)
            
            # Set z-order and alpha values for the map
            m.zorder = zorders[i]
            m.alpha = alphas[i]

            # Add map
            self._maps.append(m)

    def add_map(self, input_, zorder=None, alpha=1):
        """Adds a map to the CompositeMap
        
        Parameters
        ----------
        input_ : {sunpy.map, string}
            Map instance or filepath to map to be added
        zorder : int
            The index to use when determining where the map should lie along
            the z-axis; maps with higher z-orders appear above maps with lower
            z-orders.
        alpha : float
            Opacity at which the map should be displayed. An alpha value of 0
            results in a fully transparent image while an alpha value of 1
            results in a fully opaque image. Values between result in semi-
            transparent images.
        """
        if zorder is None:
            zorder = max([m.zorder for m in self._maps]) + 10
        
        m = BaseMap.read(input_)
        m.zorder = zorder
        m.alpha = alpha
        
        self._maps.append(m)
        
    def remove_map(self, index):
        """Removes and returns the map with the given index"""
        return self._maps.pop(index)
    
    def list_maps(self):
        """Prints a list of the currently included maps"""
        print [m.__class__ for m in self._maps]
        
    def get_alpha(self, index=None):
        """Gets the alpha-channel value for a layer in the composite image"""
        if index is None:
            return [_map.alpha for _map in self._maps]
        else:
            return self._maps[index].alpha
        
    def get_zorder(self, index = None):
        """Gets the layering preference (z-order) for a map within the
        composite.
        """
        if index is None:
            return [_map.zorder for _map in self._maps]
        else:
            return self._maps[index].zorder

    def get_colors(self, index = None):
        """Gets the colors for a map within the
        composite.
        """
        if index is None:
            return [_map.cmap for _map in self._maps]
        else:
            return self._maps[index].cmap

    def get_norm(self, index = None):
        """Gets the normalization for a map within the
        composite.
        """
        if index is None:
            return [_map.norm for _map in self._maps]
        else:
            return self._maps[index].norm

    def set_norm(self, index, norm):
        """Sets the norm for a layer in the composite image"""
        self._maps[index].norm = norm

    def set_colors(self, index, cm):
        """Sets the color map for a layer in the composite image"""
        self._maps[index].cmap = cm

    def set_alpha(self, index, alpha):
        """Sets the alpha-channel value for a layer in the composite image"""
        if 0 <= alpha <= 1:
            self._maps[index].alpha = alpha
        else:
            raise OutOfRangeAlphaValue("Alpha value must be between 0 and 1.")
        
    def set_zorder(self, index, zorder):
        """Set the layering preference (z-order) for a map within the
        composite.
        """
        self._maps[index].zorder = zorder

    def plot(self, figure=None, title="SunPy Plot", overlays=None, 
             **matplot_args):
        """Plots the composite map object using matplotlib
        
        Parameters
        ----------
        title : string
            Title to use for the plot
        overlays : list
            List of overlays to include in the plot
        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting the image.
            
        Returns
        -------
        out : matplotlib.figure.Figure
            A Matplotlib figure instance representing the composite map plot
        """
        import numpy as np

        if overlays is None:
            overlays = []

        # Create a figure and add title and axes
        if figure is None:
            figure = plt.figure()
        
        axes = figure.add_subplot(111)
        axes.set_title(title)
        
        axes.set_xlabel('X-position [' + self._maps[0].units_x + ']')
        axes.set_ylabel('Y-position [' + self._maps[0].units_y + ']')
        
        # Plot layers of composite map
        for m in self._maps:
            # Parameters for plotting
            params = {
                "origin": "lower",
                "extent": m.xrange + m.yrange,
                "cmap": m.cmap,
                "norm": m.norm(),
                "alpha": m.alpha,
                "zorder": m.zorder
            }
            params.update(matplot_args)
            
            # Use contour for contour data, and imshow otherwise
            if isinstance(m, RHESSIMap) or (m.alpha == 0):
                # Set data with values <= 0 to transparent
                contour_data = np.ma.masked_array(m, mask=(m <= 0))
                plt.contourf(contour_data, **params)
            else:
                plt.imshow(m, **params)
        
        # Adjust axes extents to include all data
        axes.axis('image')
        
        for overlay in overlays:
            figure, axes = overlay(figure, axes)

        return figure

    def show(self, title="SunPy Plot", overlays=None, **matplot_args):
        """Displays the composite map on the screen.
        
        Parameters
        ----------
        title : string
            Title to use for the plot
        overlays : list
            List of overlays to include in the plot
        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting the image.
        """
        self.plot(title, overlays, **matplot_args).show()
        
class OutOfRangeAlphaValue(ValueError):
    """Exception to raise when an alpha value outside of the range 0-1 is
    requested.
    """
    pass
