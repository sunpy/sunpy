"""
Matplotlib toolbar with extra custom actions

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

from PyQt4.QtCore import SIGNAL
from PyQt4.QtGui import QIcon, QAction
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg

class PlotToolBar(NavigationToolbar2QTAgg):
    """ Matplotlib toolbar with our own custom actions """

    def __init__(self, canvas, parent=None, coordinates=True):
        NavigationToolbar2QTAgg.__init__(self, canvas, parent, coordinates)
        self.parent = parent

    def _init_toolbar(self):
        self.actionCm_selector = self.create_action(
                        slot=self.show_cm_selector,
                        signal=SIGNAL('toggled(bool)'),
                        checkable=True,
                        icon="change_cm",
                        text="Show color map selector",
                        tip="Show color map selector"
                    )
        self.addAction(self.actionCm_selector)
        self.addSeparator()

        NavigationToolbar2QTAgg._init_toolbar(self)

        # Set status tip for all matplotlib actions
        for action in self.actions():
            action.setStatusTip(action.toolTip())
    
    def show_cm_selector(self, toggled):
        # cm_selector probably shouldn't be an attribute of a PlotToolBar instance
        # but rather of the mainwindow...
        if toggled:
            self.parent.cmComboBox.show()
            self.actionCm_selector.setToolTip("Hide color map selector")
            self.actionCm_selector.setStatusTip("Hide color map selector")
        else:
            self.parent.cmComboBox.hide()
            self.actionCm_selector.setToolTip("Show color map selector")
            self.actionCm_selector.setStatusTip("Show color map selector")

    # Move this to a gui.util module?
    def create_action(self, signal=SIGNAL('triggered()'), slot=None, text=None,
                        icon=None, tip=None, shortcut=None, checkable=None):
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
        if checkable is not None:
            action.setCheckable(checkable)

        return action

