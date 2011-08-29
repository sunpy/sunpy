from PyQt4.QtGui import QWidget, QVBoxLayout
from sunpy.gui.ui.mainwindow.widgets.figure_canvas import FigureCanvas
from sunpy.gui.ui.mainwindow.widgets.toolbars import PlotToolBar

class TabPage(QWidget):
    """ Custom widget for tab pages with canvas and plot toolbar """
    
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.canvas = FigureCanvas()
        self.plot_toolbar = PlotToolBar(self.canvas)
        layout = QVBoxLayout(self)
        layout.addWidget(self.canvas)
        layout.addWidget(self.plot_toolbar)
