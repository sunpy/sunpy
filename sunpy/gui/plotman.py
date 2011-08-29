# -*- coding:utf-8 -*-

"""
Prototype GUI Plot Manager

Plots FITS data using sunpy.Map in a Qt interface,
and provides tools for graphical plot manipulation

Notes and To-Do
===============
- Disambiguate plots with identical names (for cmselector etc.), numerical suffix (?)
- Dynamic subplot adding
- Overlaying plots, drag-drop, opacity etc.
- Restrict file types in open dialog
- Handle exceptions
- Displaying mapcube animations
- Connect to VSO and other data services
- Is gui.ui necessary?

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

import sunpy
from PyQt4.QtCore import pyqtSignature, QFileInfo
from PyQt4.QtGui import QMainWindow, QFileDialog, QVBoxLayout, QWidget
from sunpy.gui.ui.mainwindow import ui_mainwindow
from sunpy.gui.ui.mainwindow.widgets.tab_page import TabPage

class MainWindow(QMainWindow, ui_mainwindow.Ui_MainWindow):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
    
    @pyqtSignature("")
    def on_actionOpen_file_triggered(self):
        file_info = QFileInfo(QFileDialog.getOpenFileName(self, 'Open plot...'))
        file_path = str(file_info.filePath())

        tab_page = TabPage()
        tab_page.fig = sunpy.Map(file_path)
        tab_page.canvas.axes.imshow(tab_page.fig)

        tab_page.canvas.draw()
        self.tabWidget.addTab(tab_page, file_info.fileName())

if __name__ == "__main__":
    import sys
    from PyQt4.QtGui import QApplication
    app = QApplication(sys.argv)
    main = MainWindow()
    main.show()
    app.exec_()
