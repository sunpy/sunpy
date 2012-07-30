# -*- coding:utf-8 -*-
"""
SunPy PlotMan GUI

Plots FITS data using sunpy.make_map in a Qt interface,
and provides tools for graphical plot manipulation.

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""
from __future__ import absolute_import

import sunpy
from matplotlib import pyplot as plt
from PyQt4.QtCore import pyqtSignature, QFileInfo
from PyQt4.QtGui import QMainWindow, QFileDialog, QMessageBox
from sunpy.gui.ui.mainwindow import ui_mainwindow
from sunpy.gui.ui.mainwindow.widgets.tab_page import TabPage


class MainWindow(QMainWindow, ui_mainwindow.Ui_MainWindow):

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        self.colorOptionsDockWidget.hide()

    @pyqtSignature("")
    def on_actionOpen_file_triggered(self):
        file_dialog = QFileDialog(self, self.tr("Open file(s)..."),
                    filter=self.tr("FITS files (*.fit *.dst *.fits *.fts *.lilo *.lihi *.silo *.sihi *.mxlo *.mxhi *.rilo *.rihi *.vdlo *.vdhi)"))
        file_dialog.setFileMode(QFileDialog.ExistingFile)
        files = file_dialog.getOpenFileNames()
        if files:
            for file_ in files:
                file_info = QFileInfo(file_)
                file_path = str(file_info.filePath())
                self.add_tab(file_path, file_info.fileName())

    @pyqtSignature("int")
    def on_tabWidget_currentChanged(self, index):
        # index is -1 when the are no tabs
        if index >= 0:
            self.refresh_color_options()

    @pyqtSignature("int")
    def on_tabWidget_tabCloseRequested(self, index):
        self.tabWidget.removeTab(index)
        # Hide the color options dockable if last tab closed
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
    def on_cmListWidget_currentTextChanged(self, cmapname):
        self.current_tab.canvas.update_figure(cmapname=str(cmapname))

    def add_tab(self, file_path, tab_title):
        """ Adds a new tab having title 'tab_title' containing a
            TabPage widget whose FigureCanvas displays the data in 'file_path' """

        try:
            mapobject = sunpy.make_map(file_path)
            tab_page = TabPage(mapobject, self.tabWidget)
            self.tabWidget.addTab(tab_page, tab_title)
            
            # Focus new tab
            self.tabWidget.setCurrentIndex(self.tabWidget.count() - 1)
            
            # Set color options dialog appropriately
            self.initialize_color_options()
            if self.tabWidget.count() == 1:
                self.colorOptionsDockWidget.show()
        except TypeError, e:
            file_err = QMessageBox()
            file_err.setText(str(e) + '\n' + file_path)
            file_err.exec_()

    def initialize_color_options(self):
        """ Perform a first time initialisation of color
            option widgets when a new plot is opened. """
        from sunpy.cm import cm

        # Populate list widget with SunPy colormaps
        for cmap in cm.cmlist:
            self.cmListWidget.addItem(cmap)

        # Populate list widget with MPL colormaps
        self.mpl_cmaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
        for cmap in self.mpl_cmaps:
            self.cmListWidget.addItem(cmap)

        if self.current_tab.canvas.map.norm() is not None:
            # If pre-normalised, get inital clips from the matplotlib norm
            self.clipMinDoubleSpinBox.setValue(self.current_tab.canvas.map.norm().vmin)
            self.clipMaxDoubleSpinBox.setValue(self.current_tab.canvas.map.norm().vmax)
        else:
            # Otherwise, get initial clips from the map data directly.
            self.clipMinDoubleSpinBox.setValue(self.current_tab.canvas.map.min())
            self.clipMaxDoubleSpinBox.setValue(self.current_tab.canvas.map.max())

        self.scalingComboBox.setCurrentIndex(0)  # Set to linear...
        # Ideally we should set selection to the new appropriate colormap
        self.cmListWidget.clearSelection()

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
