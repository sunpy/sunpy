# -*- coding:utf-8 -*-

"""
SunPy PlotMan GUI

Plots FITS data using sunpy.Map in a Qt interface,
and provides tools for graphical plot manipulation.

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

import sunpy
from PyQt4.QtCore import pyqtSignature, QFileInfo
from PyQt4.QtGui import QMainWindow, QFileDialog
from sunpy.gui.ui.mainwindow import ui_mainwindow
from sunpy.gui.ui.mainwindow.widgets.tab_page import TabPage


class MainWindow(QMainWindow, ui_mainwindow.Ui_MainWindow):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
 
    @pyqtSignature("int")
    def on_tabWidget_tabCloseRequested(self, i):
        self.tabWidget.removeTab(i)

    @pyqtSignature("")
    def on_actionOpen_file_triggered(self):
        file_info = QFileInfo(QFileDialog.getOpenFileName(self, self.tr("Open file..."), 
                    filter=self.tr("FITS files (*.fit *.dst *.fits *.fts *.lilo *.lihi *.silo *.sihi *.mxlo *.mxhi *.rilo *.rihi *.vdlo *.vdhi)")))
        file_path = str(file_info.filePath())

        tab_page = TabPage(sunpy.Map(file_path), self.tabWidget)
        self.tabWidget.addTab(tab_page, file_info.fileName())

if __name__ == "__main__":
    import sys
    from PyQt4.QtGui import QApplication
    app = QApplication(sys.argv)
    main = MainWindow()
    main.show()
    app.exec_()
