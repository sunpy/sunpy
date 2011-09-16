#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Plot widgets for the PlotMan
subclassed from the matplotlib FigureCanvasQTAgg

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

from PyQt4.QtGui import QSizePolicy
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from matplotlib import pyplot


class FigureCanvas(FigureCanvasQTAgg):
    """ Base canvas object, resizes to fit window """

    def __init__(self, figure, parent=None):
        self.figure = figure
        FigureCanvasQTAgg.__init__(self, self.figure)
        FigureCanvasQTAgg.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)

    def update_figure(self, figure):
        pyplot.close() # Destroy any previous figures... :(
        self.figure = figure

        # This is required to get plot to size correctly...
        self.updateGeometry()
        self.resize(1, 1)

        self.draw()
