# -*- coding:utf-8 -*-

"""
SunPy PlotMan GUI

Plots FITS data using sunpy.Map in a Qt interface,
and provides tools for graphical plot manipulation.

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

import sunpy
import matplotlib.cm
from matplotlib import pyplot as plt
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
        # Focus new tab
        self.tabWidget.setCurrentIndex(self.tabWidget.count() - 1)

        self.update_color_options()

    def update_color_options(self):
        self.populate_cm_widget()
        self.set_clip_values()

    def set_clip_values(self):
        pass

    def populate_cm_widget(self):
        # Populate list widget with SunPy colormaps
        for cmap in sunpy.cm.cmlist:
            self.cmListWidget.addItem(cmap)

        # Populate list widget with MPL colormaps
        self.mpl_cmaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
        for cmap in self.mpl_cmaps:
            self.cmListWidget.addItem(cmap)
    
    @pyqtSignature("int")
    def on_gammaSlider_valueChanged(self, value):
        pass

    @pyqtSignature("QString")
    def on_cmListWidget_currentTextChanged(self, current_text):
        """ Changes color map of figure when one is selected, assumes the name
            of the colormap is identical to the name shown in the list. """

        # Check if the selected cmap belongs to MPL or SunPy and act accordingly
        cm_name = str(current_text)
        if cm_name in self.mpl_cmaps:
            new_cmap = matplotlib.cm.get_cmap(cm_name)
        else:
            new_cmap = sunpy.cm.get_cmap(cm_name)

        fig = self.current_tab._map.plot(cmap=new_cmap)
        self.current_tab.canvas.update_figure(fig)

    @property
    def current_tab(self):
        """ Return the TabPage which is currently focused """
        return self.tabWidget.currentWidget()


if __name__ == "__main__":
    import sys
    from PyQt4.QtGui import QApplication
    app = QApplication(sys.argv)
    main = MainWindow()
    main.show()
    app.exec_()
