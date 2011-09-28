"""
TabPage Widget

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

from PyQt4.QtGui import QWidget, QVBoxLayout
from sunpy.gui.ui.mainwindow.widgets.figure_canvas import FigureCanvas
from sunpy.gui.ui.mainwindow.widgets.toolbars import PlotToolBar


class TabPage(QWidget):
    """ Custom widget for tab pages with canvas and plot toolbar """

    def __init__(self, map_, parent=None):
        QWidget.__init__(self, parent)

        # Initialize canvas and toolbar
        self.canvas = FigureCanvas(map_, parent=self)
        self.plot_toolbar = PlotToolBar(self.canvas, parent=self)
       
        # Setup the page layout
        layout = QVBoxLayout(self)
        layout.addWidget(self.canvas)
        layout.addWidget(self.plot_toolbar)
