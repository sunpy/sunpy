"""A Composite Map class

Author: `Keith Hughitt <keith.hughitt@nasa.gov>`
"""
from __future__ import absolute_import

import sunpy
import numpy as np
import matplotlib.pyplot as plt

__author__ = "Keith Hughitt"
__email__ = "keith.hughitt@nasa.gov"

class CompositeMap:
    """Class representing a stack of several Maps"""    
    def __init__(self, *args):
        self._maps = []
        
        # Parse input Maps/filepaths
        zorder = 0
        
        for input_ in args:
            # Parse map
            m = sunpy.Map(input_)
            
            # Give the map a z-order and alpha
            m.zorder = zorder
            m.alpha = 1

            self._maps.append(m)
            
            zorder += 10

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
        
        m = sunpy.Map(input_)
        m.zorder = zorder
        m.alpha = alpha
        
        self._maps.append(m)
        
    def remove_map(self, index):
        """Removes and returns the map with the given index"""
        return self._maps.pop(index)
    
    def list_maps(self):
        """Prints a list of the currently included maps"""
        print [m.__class__ for m in self._maps]

    def plot(self, title="SunPy Plot", overlays=None, **matplot_args):
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
        if overlays is None:
            overlays = []

        # Create a figure and add title and axes
        fig = plt.figure()
        
        axes = fig.add_subplot(111)
        axes.set_title(title)
        
        axes.set_xlabel('X-position [' + self._maps[0].units['x'] + ']')
        axes.set_ylabel('Y-position [' + self._maps[0].units['y'] + ']')
        
        # Plot layers of composite map
        for m in self._maps:
            plt.imshow(m, origin='lower', extent=m.xrange + m.yrange, 
                       cmap=m.cmap, norm=m.norm(), alpha=m.alpha, 
                       zorder=m.zorder, **matplot_args)
        
        axes.axis('image')
        
        for overlay in overlays:
            fig, axes = overlay(fig, axes)

        return fig

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
