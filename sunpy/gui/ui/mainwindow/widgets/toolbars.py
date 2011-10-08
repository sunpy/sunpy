"""
Customizable Matplotlib toolbar

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg

class PlotToolBar(NavigationToolbar2QTAgg):

    def __init__(self, canvas, parent=None, coordinates=True):
        NavigationToolbar2QTAgg.__init__(self, canvas, parent, coordinates)
        self.parent = parent

    def _init_toolbar(self):
        NavigationToolbar2QTAgg._init_toolbar(self)

        # Set status tip for all matplotlib actions
        for action in self.actions():
            action.setStatusTip(action.toolTip())

    def home(self):
        """ Reset everything to default, including colormap when home pressed. """
        NavigationToolbar2QTAgg.home(self)
        self.parent.canvas.reset_figure() 
