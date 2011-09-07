"""
TabPage Widget

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

import matplotlib.cm
from PyQt4.QtCore import SIGNAL
from PyQt4.QtGui import QWidget, QVBoxLayout, QComboBox
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

        # Initialize color map selector
        self.cmComboBox = QComboBox()
        self.connect(self.cmComboBox,
                SIGNAL('currentIndexChanged(int)'), self.update_cm)
        self.cmComboBox.hide()
        
        # Setup the page layout
        layout = QVBoxLayout(self)
        layout.addWidget(self.canvas)
        layout.addWidget(self.plot_toolbar)
        layout.addWidget(self.cmComboBox)

    def update_cm(self, index):
        # Should this really be a method of TabPage?
        new_cmap = matplotlib.cm.get_cmap(str(self.cmComboBox.itemText(index)))
        self.canvas.figure = self._map.plot(cmap=new_cmap)
        self.canvas.draw()
