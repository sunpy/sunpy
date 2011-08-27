# -*- coding: utf-8 -*-

"""
Color Map Selector Dialog

Author: Matt Earnshaw <matt@earnshaw.org.uk>

To-do and Notes
===============
- Use sunpy.cm

"""

import matplotlib.pyplot as plt
import matplotlib.cm
from PyQt4.QtCore import pyqtSignature
from PyQt4.QtGui import QDialog
from sunpy.gui.ui.dialogs import ui_cmselector


class CMSelectorDlg(QDialog, ui_cmselector.Ui_CMSelectorDlg):

    def __init__(self, parent=None):
        QDialog.__init__(self, parent)
        self.setupUi(self)
        self.parent = parent
        self.updateUi(self.parent)

    @pyqtSignature("QString")
    def on_cmComboBox_currentIndexChanged(self, item):
        new_cmap = matplotlib.cm.get_cmap(str(item))
        self.parent.canvas.axes.imshow(self.parent.fig, cmap=new_cmap)
        self.parent.canvas.draw()

    def updateUi(self, parent):
        """ Add existing figures, subplots and color maps to
            the combo-boxes """

        # Eventually this should get a list of all canvases
        figure = self.parent.canvas         
        # Eventually this should get figure's actual name
        self.plotComboBox.addItem(str(figure)) 
        self.subplotComboBox.addItem(str(figure.axes.get_axes()))
        maps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
        for map in maps:
            self.cmComboBox.addItem(map)

if __name__ == "__main__":
    import sys
    from PyQt4.QtGui import QApplication
    app = QApplication(sys.argv)
    dlg = CMSelectorDlg()
    dlg.show()
    app.exec_()
