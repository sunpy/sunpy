#-*- coding: utf-8 -*-
# Author: Matt Earnshaw <matt@earnshaw.org.uk>

from __future__ import absolute_import
import os
import sys
import sunpy
from PyQt4.QtGui import QApplication
from sunpy.gui.mainwindow import MainWindow
from sunpy.io import UnrecognizedFileTypeError 

class Plotman(object):
    """ Wraps a MainWindow so PlotMan instances can be created via the CLI.

        Example
        -------
        from sunpy.gui import Plotman
        plots = Plotman("data/examples")
        plots.show()
    """

    def __init__(self, *paths):
        """ *paths: directories containing FITS paths
                   or FITS paths to be opened in PlotMan """
        self.app = QApplication(sys.argv)
        self.main = MainWindow()
        self.open_paths(paths)

    def open_paths(self, paths):
        for path in paths:
            if os.path.isfile(path):
                try:
                    self.main.add_tab(path, os.path.basename(path))
                except UnrecognizedFileTypeError:
                    pass
            elif os.path.isdir(path):
                for file_ in os.listdir(path):
                    if path[-1] != os.path.sep:
                        path += os.path.sep
                    try:
                        self.main.add_tab(path + file_, file_)
                    except UnrecognizedFileTypeError:
                        pass
            else:
                raise IOError("Path " + path + " does not exist.")

    def show(self):
        self.main.show()
        self.app.exec_()

if __name__=="__main__":
    from sunpy.gui import Plotman
    plots = Plotman(sunpy.AIA_171_IMAGE)
    plots.show()
