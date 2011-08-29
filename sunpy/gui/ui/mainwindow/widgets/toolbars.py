from PyQt4.QtCore import SIGNAL
from PyQt4.QtGui import QIcon, QAction
from sunpy.gui.ui.dialogs import cmselector
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg

class PlotToolBar(NavigationToolbar2QTAgg):
    """ Matplotlib toolbar with our own custom actions """

    def __init__(self, canvas, parent=None, coordinates=True):
        NavigationToolbar2QTAgg.__init__(self, canvas, parent, coordinates)

    def _init_toolbar(self):
        NavigationToolbar2QTAgg._init_toolbar(self)
        change_cm = self.create_action(
                        slot=self.change_cm,
                        icon="change_cm",
                        tip="Change plot colour map..."
                    )
        self.addAction(change_cm)
    
    def change_cm(self):
        # cm_selector probably shouldn't be an attribute of a PlotToolBar instance
        # but rather of the mainwindow...
        if not hasattr(self, "cm_selector"):
            parent = self.window()
            self.cm_selector = cmselector.CMSelectorDlg(parent)
        self.cm_selector.show()
        self.cm_selector.raise_()
        self.cm_selector.activateWindow()
    
    def create_action(self, signal=SIGNAL('triggered()'), slot=None, text=None,
                        icon=None, tip=None, shortcut=None):
        """ Helper function for creating useful QAction objects"""

        action = QAction(self)
        if slot is not None:
            self.connect(action, signal, slot)
        if text is not None:
            action.setText(text)
        if tip is not None:
            action.setStatusTip(tip)
            action.setToolTip(tip)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)

        return action

