#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Plot widgets for the PlotMan
subclassed from the matplotlib FigureCanvasQTAgg

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

import sunpy
import matplotlib
from PyQt4.QtGui import QSizePolicy
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg

class FigureCanvas(FigureCanvasQTAgg):
    """ General plot widget, resizes to fit window """

    def __init__(self, map_, parent=None):
        self.parent = parent
        self.map_ = map_
        self.figure = self.map_.plot()

        # Code duplication from plotman.py!
        if self.map_.norm() is not None:
            # If pre-normalised, get inital clips from the matplotlib norm
            self.vmin = self.map_.norm().vmin
            self.vmax = self.map_.norm().vmax
        else:
            # Otherwise, get initial clips from the map data directly.
            self.vmin = self.map_.min()
            self.vmax = self.map_.max()
        
        self.scaling = "Linear" # Shouldn't be a default.

        # Matplotlib kwargs
        self.params = {
                        "cmap": self.map_.cmap,
                        "norm": self.map_.norm(),
                      }

        FigureCanvasQTAgg.__init__(self, self.figure)
        FigureCanvasQTAgg.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

    def get_cmap(self, cmap_name, gamma=None):
        """ Given a sunpy or matplotlib cmap name, returns the cmap. """
        mpl_cmaps = sorted(m for m in matplotlib.pyplot.cm.datad
                            if not m.endswith("_r"))
        if cmap_name in mpl_cmaps:
            return matplotlib.cm.get_cmap(cmap_name)
        else:
            return sunpy.cm.get_cmap(cmap_name)

    def get_norm(self):
        """ Return a matplotlib Normalize according to clip values and scale type """
        if self.vmin < self.vmax:
            if self.scaling == "Linear": # This is not robust! 
                return matplotlib.colors.Normalize(vmin=self.vmin,
                                                   vmax=self.vmax)
            elif self.scaling == "Logarithmic":
                if self.vmin <= 0:
                    # Cannot log scale 0 or negative data
                    self.window().colorErrorLabel.setText("Cannot log scale zero or negative data!")
                    return self.params["norm"]
                else:
                    return matplotlib.colors.LogNorm(vmin=self.vmin,
                                                     vmax=self.vmax)
        else:
            self.window().colorErrorLabel.setText("Min clip value must not exceed max clip value!")
            return self.params["norm"]


    def update_figure(self, cmap_name=None, vmin=None, vmax=None, scaling=None):
        # Destroy any previous figures
        matplotlib.pyplot.close()

        # Clear error messages
        self.window().colorErrorLabel.clear()        

        if cmap_name is not None:
            self.params.update({"cmap": self.get_cmap(cmap_name)})
        elif vmax is not None:
            self.vmax = vmax
            self.params.update({"norm": self.get_norm()})
        elif vmin is not None:
            self.vmin = vmin
            self.params.update({"norm": self.get_norm()})
        elif scaling is not None:
            self.scaling = scaling
            self.params.update({"norm": self.get_norm()})

        self.figure = self.map_.plot(**self.params)
        self.resize_figure() 

    def reset_figure(self):
        self.figure = self.map_.plot()
        self.resize_figure()

    def resize_figure(self):
        self.updateGeometry()
        self.resize(1, 1) # It works, OK?
        self.draw()

    @property
    def cmap(self):
        return self.params["cmap"]

    @property
    def norm(self):
        return self.params["norm"]
