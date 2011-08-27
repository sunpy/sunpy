"""
Prototype GUI Plot Manager

Plots FITS data using sunpy.Map in a Qt interface,
and provides tools for graphical plot manipulation

Notes and To-Do
===============
- Tabbed interface for seperate FigureCanvases (?)
- Dynamic subplot adding
- Overlaying plots, drag-drop, opacity etc.
- Restrict file types in open dialog
- Handle exceptions
- Colour-map selector window/interface with sunpy.cm
- Displaying mapcube animations
- Connect to VSO and other data services

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

from PyQt4.QtCore import pyqtSignature
from PyQt4.QtGui import QMainWindow, QFileDialog
from sunpy.gui.ui.mainwindow import ui_mainwindow
from sunpy.gui.dialogs import cmselector

class MainWindow(QMainWindow, ui_mainwindow.Ui_MainWindow):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)

    @pyqtSignature("")
    def on_actionOpen_file_triggered(self):
        file_path = str(QFileDialog.getOpenFileName(self, 'Open plot...'))

    @pyqtSignature("")
    def on_actionChange_color_map_triggered(self):
        cm_selector = cmselector.CMSelectorDlg(self)
        cm_selector.exec_()

if __name__ == "__main__":
    import sys
    from PyQt4.QtGui import QApplication
    app = QApplication(sys.argv)
    main = MainWindow()
    main.show()
    app.exec_()
