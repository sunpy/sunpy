#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Plot widgets for the PlotMan
subclassed from the matplotlib FigureCanvasQTAgg

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

from PyQt4.QtGui import QSizePolicy
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg


class FigureCanvas(FigureCanvasQTAgg):
    """ Base canvas object, resizes to fit window """

    def __init__(self, figure, parent=None):
        self.figure = figure
        FigureCanvasQTAgg.__init__(self, self.figure)
        FigureCanvasQTAgg.setSizePolicy(self,
                                   QSizePolicy.MinimumExpanding,
                                   QSizePolicy.MinimumExpanding)
        FigureCanvasQTAgg.updateGeometry(self)

    def update_figure(self, figure):
        self.figure = figure
        self.updateGeometry()
        self.resize(1, 1)
        self.draw()
