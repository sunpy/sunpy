# -*- coding: utf-8 -*-

"""
Color Map Selector Dialog

Author: Matt Earnshaw <matt@earnshaw.org.uk>

To-do and Notes
===============
- Clean up horrible mix of mainWindow, parent etc.
- Generally make more robust and less verbose
- Disambiguate plots with identical names (probably handle this in mainwindow)
- Remove plot from ui. Only needs to refer to subplots.
- Add "all" option to subplots list
- Use sunpy.cm (?)
- Add subplots functionality

"""

import matplotlib.pyplot as plt
import matplotlib.cm
from PyQt4.QtCore import pyqtSignature, QString
from PyQt4.QtGui import QDialog
from sunpy.gui.ui.dialogs import ui_cmselector


class CMSelectorDlg(QDialog, ui_cmselector.Ui_CMSelectorDlg):

    def __init__(self, parent=None):
        QDialog.__init__(self, parent)
        self.setupUi(self)
        self.mainWindow = parent
        self.updateUi(parent)

    @pyqtSignature("int")
    def on_cmComboBox_currentIndexChanged(self, i):
        # Note: assumes that item at combobox index corresponds to tab_page at same index.
        tab_page = self.mainWindow.tabWidget.widget(self.plotComboBox.currentIndex())
        new_cmap = matplotlib.cm.get_cmap(str(self.cmComboBox.itemText(i)))
        tab_page.canvas.axes.imshow(tab_page.fig, cmap=new_cmap)
        tab_page.canvas.draw()

    def updateUi(self, parent):
        """ Add existing subplots and color maps to the combo-boxes """

        # Assumes that each tab = canvas
        tabWidget = parent.tabWidget
        for i in range(tabWidget.count()):
            self.plotComboBox.addItem(tabWidget.tabText(i))
            # QString(str()) ... ugly.
            self.subplotComboBox.addItem(QString(str(tabWidget.widget(i).canvas.axes.get_axes())))

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
