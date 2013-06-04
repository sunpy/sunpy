#-*- coding: utf-8 -*-
# Author: Matt Earnshaw <matt@earnshaw.org.uk>

from __future__ import absolute_import
import os
import sys
import sunpy
from PyQt4.QtGui import QApplication
from sunpy.gui.mainwindow import MainWindow
from sunpy.io.file_tools import UnrecognizedFileTypeError

class Plotman(object):
    """ Wraps a MainWindow so PlotMan instances can be created via the CLI.

        Examples
        --------
        from sunpy.gui import Plotman
        plots = Plotman("data/examples")
        plots.show()
    """

    def __init__(self, *paths):
        """ *paths: directories containing FITS paths
                   or FITS paths to be opened in PlotMan """
        self.app = QApplication(sys.argv)
        self.main = MainWindow()
        self.open_files(paths)

    def open_files(self, inputs):
        VALID_EXTENSIONS = [".jp2", ".fits", ".fts"]
        
        to_open = []
        
        # Determine files to process
        for input_ in inputs:
            if os.path.isfile(input_):
                to_open.append(input_)
            elif os.path.isdir(input_):
                for file_ in os.listdir(input_):
                    to_open.append(file_)
            else:
                raise IOError("Path " + input_ + " does not exist.")

        # Load files
        for filepath in to_open:
            name, ext = os.path.splitext(filepath) #pylint: disable=W0612
            
            if ext.lower() in VALID_EXTENSIONS:
                try:
                    self.main.add_tab(filepath, os.path.basename(filepath))
                except UnrecognizedFileTypeError:
                    pass

    def show(self):
        self.main.show()
        self.app.exec_()

if __name__=="__main__":
    from sunpy.gui import Plotman
    plots = Plotman(sunpy.AIA_171_IMAGE)
    plots.show()
