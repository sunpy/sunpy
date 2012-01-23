#!/usr/bin/env python
#-*- coding:utf-8 -*-
import sys
import datetime
import sunpy
import Ui_RGBComposite
from sunpy.net import helioviewer as hv
from PyQt4 import QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

def main():
    app = QtGui.QApplication(sys.argv)
    win = RGBCompositeImageApp()
    win.show()
    sys.exit(app.exec_())

#        super(RGBCompositeImageApp, self).__init__()
class RGBCompositeImageApp(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.ui = Ui_RGBComposite.Ui_RGBComposite()
        self.ui.setupUi(self)
        
        # Loaded images
        self.red = None
        self.green = None
        self.blue = None
        
        # Setup UI
        self._load_data_sources()
        self._load()
        
    def _load_data_sources(self):
        """Downloads and displays latest images for default wavelengths"""
        self._datasources = hv.get_data_sources()['SDO']['AIA']['AIA']
        sorted_datasources = sorted(self._datasources, key=int)
        
        for wl in sorted_datasources:
            self.ui.redWavelengthSelect.addItem(wl, self._datasources[wl])
            self.ui.greenWavelengthSelect.addItem(wl, self._datasources[wl])
            self.ui.blueWavelengthSelect.addItem(wl, self._datasources[wl])
            
        # Default wavelengths: 304, 193, 171
        self.ui.redWavelengthSelect.setCurrentIndex(5)
        self.ui.redWavelengthSelect.setCurrentIndex(3)
        self.ui.redWavelengthSelect.setCurrentIndex(2)
        
    def _load(self):
        """Finish UI initialization and load default RGB composite image"""
        # Default to the current date
        self.ui.dateTimeEdit.setDateTime(datetime.datetime.utcnow())
        
        now = datetime.datetime.utcnow()
        
        # Load initial RGB image
        self.red = sunpy.make_map(hv.get_jp2_image(now, sourceId=self._datasources['304']['sourceId']))
        self.green = sunpy.make_map(hv.get_jp2_image(now, sourceId=self._datasources['193']['sourceId']))
        self.blue = sunpy.make_map(hv.get_jp2_image(now, sourceId=self._datasources['171']['sourceId']))
        
        self.ui.redPreviewImage = Miniplot(self.red)
        
class Miniplot(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, map_, parent=None, width=1.28, height=1.28, dpi=100):
        self.figure = map_.plot()
        FigureCanvas.__init__(self, self.figure)
        FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

if __name__ == '__main__':
    sys.exit(main())