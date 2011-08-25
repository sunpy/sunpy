#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Plot widgets for the PlotMan
Subclassed from the matplotlib FigureCanvasQTAgg 
Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

from matplotlib.figure import Figure
from PyQt4.QtGui import QSizePolicy
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg

class BasePlot(FigureCanvasQTAgg):
    """ Base plot object, resizes to fit window """

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.axes.hold(False)
        FigureCanvasQTAgg.__init__(self, self.fig)
        FigureCanvasQTAgg.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvasQTAgg.updateGeometry(self)



