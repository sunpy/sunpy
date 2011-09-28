# -*- coding:utf-8 -*-

from __future__ import absolute_import

"""
SunPy PlotMan GUI

Plots FITS data using sunpy.Map in a Qt interface,
and provides tools for graphical plot manipulation.

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

import sunpy
from matplotlib import pyplot as plt
from PyQt4.QtCore import pyqtSignature, QFileInfo
from PyQt4.QtGui import QMainWindow, QFileDialog
from sunpy.gui.ui.mainwindow import ui_mainwindow
from sunpy.gui.ui.mainwindow.widgets.tab_page import TabPage


class MainWindow(QMainWindow, ui_mainwindow.Ui_MainWindow):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.colorOptionsDockWidget.hide()

    @pyqtSignature("")
    def on_actionOpen_file_triggered(self):
        file_info = QFileInfo(QFileDialog.getOpenFileName(self, self.tr("Open file..."), 
                    filter=self.tr("FITS files (*.fit *.dst *.fits *.fts *.lilo *.lihi *.silo *.sihi *.mxlo *.mxhi *.rilo *.rihi *.vdlo *.vdhi)")))
        if file_info.fileName():
            file_path = str(file_info.filePath())
            self.add_tab(sunpy.Map(file_path), file_info.fileName())

    @pyqtSignature("int")
    def on_tabWidget_currentChanged(self):
        self.refresh_color_options()

    @pyqtSignature("int")
    def on_tabWidget_tabCloseRequested(self, i):
        self.tabWidget.removeTab(i)
        if self.tabWidget.count() == 0:
            self.colorOptionsDockWidget.hide()

    @pyqtSignature("")
    def on_clipMinDoubleSpinBox_editingFinished(self):
        vmin = self.clipMinDoubleSpinBox.value()
        self.current_tab.canvas.update_figure(vmin=vmin)
        self.refresh_color_options()

    @pyqtSignature("")
    def on_clipMaxDoubleSpinBox_editingFinished(self):
        vmax = self.clipMaxDoubleSpinBox.value()
        self.current_tab.canvas.update_figure(vmax=vmax)
        self.refresh_color_options()

    @pyqtSignature("QString")
    def on_scalingComboBox_currentIndexChanged(self, scaling):
        self.current_tab.canvas.update_figure(scaling=scaling)
        self.refresh_color_options()
  
    @pyqtSignature("QString")
    def on_cmListWidget_currentTextChanged(self, cmap_name):
        self.current_tab.canvas.update_figure(cmap_name=str(cmap_name))

    def add_tab(self, map_object, tab_title):
        """ Adds a new tab having title 'tab_title' containing a
            TabPage widget whose FigureCanvas displays 'map_object' """
        tab_page = TabPage(map_object, self.tabWidget)
        self.tabWidget.addTab(tab_page, tab_title)
        # Focus new tab
        self.tabWidget.setCurrentIndex(self.tabWidget.count() - 1)
        # Set color options dialog appropriately
        self.initialize_color_options()
        if self.tabWidget.count() == 1:
            self.colorOptionsDockWidget.show()
          
    def initialize_color_options(self):
        """ Perform a first time initialisation of color 
            option widgets when a new plot is opened. """
        
        # Populate list widget with SunPy colormaps
        for cmap in sunpy.cm.cmlist:
            self.cmListWidget.addItem(cmap)

        # Populate list widget with MPL colormaps
        self.mpl_cmaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
        for cmap in self.mpl_cmaps:
            self.cmListWidget.addItem(cmap)

        if self.current_tab.canvas.map_.norm is not None:
            # If pre-normalised, get inital clips from the matplotlib norm
            self.clipMinDoubleSpinBox.setValue(self.current_tab.canvas.map_.norm.vmin)
            self.clipMaxDoubleSpinBox.setValue(self.current_tab.canvas.map_.norm.vmax)
        else:
            # Otherwise, get initial clips from the map data directly.
            self.clipMinDoubleSpinBox.setValue(self.current_tab.canvas.map_.min())
            self.clipMaxDoubleSpinBox.setValue(self.current_tab.canvas.map_.max())

    # This method should be combined with the above to avoid code duplication.
    def refresh_color_options(self):
        """ Set widgets according to focused plot's properties. """
        self.clipMinDoubleSpinBox.setValue(self.current_tab.canvas.vmin)
        self.clipMaxDoubleSpinBox.setValue(self.current_tab.canvas.vmax)

        # ... NOT a good solution
        if self.current_tab.canvas.scaling == "Linear":
            self.scalingComboBox.setCurrentIndex(0)
        else:
            self.scalingComboBox.setCurrentIndex(1)
 
    @property
    def current_tab(self):
        return self.tabWidget.currentWidget()


if __name__ == "__main__":
    import sys
    from PyQt4.QtGui import QApplication
    app = QApplication(sys.argv)
    main = MainWindow()
    main.show()
    app.exec_()
