"""
TabPage Widget

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

import matplotlib.pyplot as plt
import matplotlib.cm
from PyQt4.QtCore import SIGNAL
from PyQt4.QtGui import QWidget, QVBoxLayout, QComboBox
from sunpy.gui.ui.mainwindow.widgets.figure_canvas import FigureCanvas
from sunpy.gui.ui.mainwindow.widgets.toolbars import PlotToolBar

class TabPage(QWidget):
    """ Custom widget for tab pages with canvas and plot toolbar """
    
    def __init__(self, figure, parent=None):
        QWidget.__init__(self, parent)
        self.fig = figure
        self.canvas = FigureCanvas(parent=self)
        self.plot_toolbar = PlotToolBar(self.canvas, parent=self)
        # color map selector
        self.cmComboBox = QComboBox(parent=self)
        self.connect(self.cmComboBox, SIGNAL('currentIndexChanged(int)'), self.update_cm)
        self.cmComboBox.hide()
        layout = QVBoxLayout(self)
        layout.addWidget(self.canvas)
        layout.addWidget(self.plot_toolbar)
        layout.addWidget(self.cmComboBox)
        self.updateUi()

    def update_cm(self, index):
        # Should this be a method of TabPage?
        new_cmap = matplotlib.cm.get_cmap(str(self.cmComboBox.itemText(index)))
        self.canvas.axes.imshow(self.fig, cmap=new_cmap)
        self.canvas.draw()

    def updateUi(self):
        # Get all mpl color maps and add to comboBox
        maps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
        for map in maps:
            self.cmComboBox.addItem(map)
