#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Pathfinder for a prototype plot manager for SunPy
("A first program in PyQt4")

Plots FITS data using sunpy.Map in a Qt interface,
and allows the user to change the color-map of the plot

Notes and To-Do
===============
- Get addActions working
- Use the builtin matplotlib navigation toolbar?

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

import os
import sys
import sunpy
from PyQt4.QtGui import QApplication, QMainWindow, QAction, QIcon, QFileDialog
from PyQt4.QtCore import SIGNAL, SLOT
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class MainWindow(QMainWindow):

    def __init__(self, parent=None):
        """ Initialise the application's main window and its widgets. """

        QMainWindow.__init__(self, parent)
        self.setGeometry(0, 0, 800, 800)
        self.setWindowTitle("SunPy PlotMan GUI - Pathfinding Prototype")

        status_bar = self.statusBar()
        status_bar.showMessage("Ready.")

        # Setup File actions toolbar
        file_toolbar = self.addToolBar("File")
        open_plot = self.create_action(slot=self.open_plot, icon_path="icons/open_plot.png",
                                       tip="Open a FITS file for plotting...")
        save_plot = self.create_action(slot=self.save_plot, icon_path="icons/save_plot.png",
                                       tip="Save plot to PNG file...")
        exit_app = self.create_action(slot=SLOT('close()'), icon_path="icons/exit.png",
                                      tip="Exit SunPy PlotMan GUI.")
        # Should be able to use addActions here... couldn't get it to work ATM.
        file_toolbar.addAction(open_plot)
        file_toolbar.addAction(save_plot)
        file_toolbar.addAction(exit_app)

        # Setup Plot actions toolbar
        plot_toolbar = self.addToolBar("Plot")
        change_cm = self.create_action(slot=self.change_cm, icon_path="icons/change_cm.png",
                    tip="Change the Map's color-map")
        plot_toolbar.addAction(change_cm)

    def create_action(self, signal=SIGNAL('triggered()'), slot=None, text=None,
                        icon_path=None, tip=None, shortcut=None):
        """ Helper function for creating useful QAction objects"""

        action = QAction(self)
        if slot is not None:
            self.connect(action, signal, slot)
        if text is not None:
            action.setText(text)
        if tip is not None:
            action.setStatusTip(tip)
            action.setToolTip(tip)
        if icon_path is not None:
            if os.path.exists(icon_path):
                action.setIcon(QIcon(icon_path))
            else:
                IOError("Icon could not be set, " + icon_path + " does not exist")
        if shortcut is not None:
            action.setShortcut(shortcut)
        return action

    def open_plot(self):
        """docstring for open_plot"""
        file_path = QFileDialog.getOpenFileName(self, 'Open plot...', '/home')
        figure = sunpy.Map(file_path).plot(show_plot=False)
        self.canvas = FigureCanvas(figure)
        self.setCentralWidget(self.canvas)

    def save_plot(self):
        """docstring for save_plot"""
        pass

    def change_cm(self):
        """docstring for change_cm"""
        pass

    def update_plot(self):
        """docstring for update_plot"""
        pass

if __name__ == "__main__":
    app = QApplication(sys.argv)
    main = MainWindow()
    main.show()
    app.exec_()
