#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Pathfinder for a prototype plot manager for SunPy
("A first program in PyQt4")

Plots FITS data using sunpy.Map in a Qt interface,
and allows the user to change the color-map of the plot

Notes and To-Do
===============
- Add described functionality
- Get addActions working
- Use the builtin matplotlib navigation toolbar?
- More file formats for saving
- Restrict file types in open dialog
- Handle expections

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

import os
import sys

sys.path.append("../../..") # urgh
import sunpy
import qrc_resources
from PyQt4.QtGui import *
from PyQt4.QtCore import SIGNAL, SLOT
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.cm as cm


class PlotWidget(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.axes.hold(False)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


class MainWindow(QMainWindow):

    def __init__(self, parent=None):
        """ Initialise the application's main window and its widgets. """

        QMainWindow.__init__(self, parent)
        self.setGeometry(0, 0, 800, 800)
        self.setWindowTitle("SunPy PlotMan GUI - Pathfinding Prototype -")

        # Setup matplotlib figure canvas with test data
        self.canvas = PlotWidget()
        self.setCentralWidget(self.canvas)

        # Setup status bar
        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready.")

        # Setup File actions toolbar
        file_toolbar = self.addToolBar("File")
        file_toolbar_actions = QActionGroup(self)
        open_plot = self.create_action(slot=self.open_plot, icon="open_plot",
                                       tip="Open a FITS file for plotting...")
        save_plot = self.create_action(slot=self.save_plot, icon="save_plot",
                                       tip="Save plot to PNG file...")
        exit_app = self.create_action(slot=SLOT('close()'), icon="exit",
                                      tip="Exit SunPy PlotMan GUI.")
        # Should be able to use addActions here... couldn't get it to work ATM.
        file_toolbar.addAction(open_plot)
        file_toolbar.addAction(save_plot)
        file_toolbar.addAction(exit_app)

        # Setup Plot actions toolbar
        plot_toolbar = self.addToolBar("Plot")
        change_cm = self.create_action(slot=self.change_cm, icon="change_cm",
                                       tip="Change the Map's color-map")
        plot_toolbar.addAction(change_cm)

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
        """docstring for open_plot"""
        file_path = str(QFileDialog.getOpenFileName(self, 'Open plot...'))
        self.fig = sunpy.Map(file_path)
        self.canvas.axes.imshow(self.fig)
        self.canvas.draw()

    def save_plot(self):
        """docstring for save_plot"""
        file_path = str(QFileDialog.getSaveFileName(self, 'Save plot...'))
        self.canvas.fig.savefig(file_path, format="png")
        self.status_bar.showMessage("Saved.", 3000)

    def change_cm(self):
        """docstring for change_cm"""
        self.canvas.axes.imshow(self.fig, cmap=cm.hot)
        self.canvas.draw()
        self.status_bar.showMessage("Color-map changed.", 3000)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    main = MainWindow()
    main.show()
    app.exec_()
