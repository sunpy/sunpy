#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Prototype GUI Plot Manager

Plots FITS data using sunpy.Map in a Qt interface,
and provides tools for graphical plot manipulation

Notes and To-Do
===============
- Migrate completely to using Qt Designer
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

import sys
import sunpy
from resources import qrc_resources
from PyQt4.QtGui import (QApplication, QWidget, QMainWindow,
                        QAction, QIcon, QFileDialog, QVBoxLayout)
from PyQt4.QtCore import SIGNAL, SLOT, QSettings, QVariant
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg
import matplotlib.cm as cm
from widgets.figure_canvas import BaseFigureCanvas
from dialogs import cmselector

class MainWindow(QMainWindow):

    def __init__(self, parent=None):
        """ Initialise the application's main window and its widgets. """

        QMainWindow.__init__(self, parent)
        settings = QSettings()
        self.restoreGeometry(settings.value("Geometry").toByteArray())
        self.setWindowTitle(QApplication.applicationName())

        # Setup central widget and layout
        self.main_widget = QWidget(self)
        layout = QVBoxLayout(self.main_widget)
        self.setCentralWidget(self.main_widget)
        self.canvas = BaseFigureCanvas()
        layout.addWidget(self.canvas)

        # Setup status bar
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready.")

        # Setup toolbars
        # Instantiate a MPL built-in toolbar so we can reuse the actions
        mpl_toolbar = NavigationToolbar2QTAgg(self.canvas, None)
        self.general_toolbar = self.addToolBar("General")
        self.plot_toolbar = self.addToolBar("Plot")

        open_plot = self.create_action(
                        slot=self.open_plot,
                        icon="open_plot",
                        tip="Open a FITS file for plotting..."
                    )
        save_plot = self.create_action(
                        slot=mpl_toolbar.save_figure,
                        icon="save_plot",
                        tip="Save plot to PNG..."
                    )
        home = self.create_action(
                        slot=mpl_toolbar.home,
                        icon="home",
                        tip="Reset original view."
                    )
        forward = self.create_action(
                        slot=mpl_toolbar.forward,
                        icon="forward",
                        tip="Forward to next view."
                    )
        back = self.create_action(
                        slot=mpl_toolbar.back,
                        icon="back",
                        tip="Back to previous view."
                    )
        pan = self.create_action(
                        slot=mpl_toolbar.pan,
                        icon="move",
                        tip="Pan axes with left mouse, zoom with right."
                    )
        config_subplots = self.create_action(
                        slot=mpl_toolbar.configure_subplots,
                        icon="subplots",
                        tip="Configure subplots."
                    )
        zoom = self.create_action(
                        slot=mpl_toolbar.zoom,
                        icon="zoom_to_rect",
                        tip="Zoom to rectangle."
                    )
        change_cm = self.create_action(
                        slot=self.change_cm,
                        icon="change_cm",
                        tip="Change plot colour map..."
                    )
        exit = self.create_action(
                        slot=SLOT('close()'),
                        icon="exit",
                        tip="Exit SunPy PlotMan."
                    )

        self.general_toolbar.addActions([open_plot, exit])
        self.plot_toolbar.addActions([home, back, forward, pan, config_subplots, zoom, change_cm, save_plot])
        self.updateUi()

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

    def open_plot(self):
        """ Open a file to plot """
        file_path = str(QFileDialog.getOpenFileName(self, 'Open plot...'))
        self.fig = sunpy.Map(file_path)
        self.canvas.axes.imshow(self.fig)
        self.canvas.draw()
        self.updateUi()

    def change_cm(self):
        """ Display color map selector dialog """
        cm_selector = cmselector.CMSelectorDlg(self)
        cm_selector.exec_()

    def updateUi(self):
        self.plot_toolbar.setEnabled(hasattr(self, "fig"))

    def closeEvent(self, event):
        """ Save the current settings (window size, position, etc.) when the GUI is closed. """
        settings = QSettings()
        settings.setValue("Geometry", QVariant(self.saveGeometry()))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setOrganizationName("SunPy")
    app.setOrganizationDomain("sunpy.org")
    app.setApplicationName("SunPy Plot Manager")
    app.setWindowIcon(QIcon(":/main.png"))
    main = MainWindow()
    main.show()
    app.exec_()
