"""
TabPage Widget

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

from PyQt4.QtGui import QWidget, QVBoxLayout
from sunpy.gui.ui.mainwindow.widgets.figure_canvas import FigureCanvas
from sunpy.gui.ui.mainwindow.widgets.toolbars import PlotToolBar


class TabPage(QWidget):
    """ Custom widget for tab pages with canvas and plot toolbar """

    def __init__(self, _map, parent=None):
        """ 
        _map: sunpy.Map object
        """
        QWidget.__init__(self, parent)

        # Initialize canvas and toolbar
        self._map = _map
        self.figure = _map.plot()
        self.canvas = FigureCanvas(self.figure)

        # We pass self to PlotToolBar so it becomes a parent
        # of TabPage and can refer to its widgets... clumsy?
        self.plot_toolbar = PlotToolBar(self.canvas, self)
       
        # Setup the page layout
        layout = QVBoxLayout(self)
        layout.addWidget(self.canvas)
        layout.addWidget(self.plot_toolbar)
