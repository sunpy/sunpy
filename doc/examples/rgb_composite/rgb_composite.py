#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
RGB Composite Image Demo

Last Update: Jan 25, 2012
keith.hughitt@nasa.gov

This example application demonstrates how SunPy can be used to build GUI
applications using PyQt. The purpose of this simple application is to
combine three JPEG 2000 images from Helioviewer.org into a single composite
RGB image and provide a mechanism to control the contribution of each color
channel to the final image.

The GUI was built using Qt 4 Designer. To generate the Python code associated
with rgb_composite.ui, use the pyuic4 tool, e.g.:

  pyuic4 rgb_composite.ui > Ui_RGBComposite.py
  
TODO:
    * Fix bug causing composite image plot to become distored when channel
      weights are adjusted.
    * Have file -> save call savefig() on the rgb image
    * Show actual dates below each image.
    * Wavelength/date adjustment support.
    * Refactor/simplify
    * Make all images expand to fill available space.
"""
import sys
import datetime
import sunpy
import Ui_RGBComposite
import numpy as np
from sunpy.net import helioviewer as hv
from sunpy.map import BaseMap
from sunpy.util.util import toggle_pylab
from PyQt4 import QtGui, QtCore
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

def main():
    app = QtGui.QApplication(sys.argv)
    win = RGBCompositeImageApp()
    win.show()
    sys.exit(app.exec_())

class RGBCompositeImageApp(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.ui = Ui_RGBComposite.Ui_RGBComposite()
        self.ui.setupUi(self)
        
        # Loaded images
        self.red = None
        self.green = None
        self.blue = None
        
        # Color channel weights
        self._weights=[1., 1., 1.]
        
        # Setup UI
        self._load_data_sources()
        self._load()
        
        # Initialize event handlers
        self._initEvents()
        
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
        self.ui.greenWavelengthSelect.setCurrentIndex(3)
        self.ui.blueWavelengthSelect.setCurrentIndex(2)
        
    def _load(self):
        """Finish UI initialization and load default RGB composite image"""
        # Default to the current date
        self.ui.dateTimeEdit.setDateTime(datetime.datetime.utcnow())
        
        now = datetime.datetime.utcnow()
        
        # Load initial RGB image
        self.red = sunpy.make_map(hv.get_jp2_image(now, sourceId=self._datasources['304']['sourceId']))
        self.green = sunpy.make_map(hv.get_jp2_image(now, sourceId=self._datasources['193']['sourceId']))
        self.blue = sunpy.make_map(hv.get_jp2_image(now, sourceId=self._datasources['171']['sourceId']))
        
        self._updateRedPreview()
        self._updateGreenPreview()
        self._updateBluePreview()
        self._updateCompositeImage()
        
    def _updateCompositeImage(self):
        """Updates the RGB composite image"""
        rgb = RGBCompositeMap(self.red, self.green, self.blue)
        self.ui.compositeContainer.removeWidget(self.ui.compositePlaceholder)
        self.ui.compositePlaceholder.close()
        self.ui.compositeImage = SunPyPlot(rgb, 512, 512)
        self.ui.compositeContainer.addWidget(self.ui.compositeImage, 1)
        self.ui.compositeContainer.update()

    def _updateRedPreview(self):
        """Updates the red preview image"""
        self.ui.redPreview.removeWidget(self.ui.redPlaceholder)
        self.ui.redPlaceholder.close()
        self.ui.redPreviewImage = SunPyPlot(self.red, 256, 256, cmap=cm.Reds_r)
        self.ui.redPreview.addWidget(self.ui.redPreviewImage, 1)
        self.ui.redPreview.update()
        
    def _updateGreenPreview(self):
        """Updates the green preview image"""
        self.ui.greenPreview.removeWidget(self.ui.greenPlaceholder)
        self.ui.greenPlaceholder.close()
        self.ui.greenPreviewImage = SunPyPlot(self.green, 256, 256, cmap=cm.Greens_r)
        self.ui.greenPreview.addWidget(self.ui.greenPreviewImage, 1)
        self.ui.greenPreview.update()
        
    def _updateBluePreview(self):
        """Updates the blue preview image"""
        self.ui.bluePreview.removeWidget(self.ui.bluePlaceholder)
        self.ui.bluePlaceholder.close()
        self.ui.bluePreviewImage = SunPyPlot(self.blue, 256, 256, cmap=cm.Blues_r)
        self.ui.bluePreview.addWidget(self.ui.bluePreviewImage, 1)
        self.ui.bluePreview.update()
        
    def _initEvents(self):
        """Initialize event handlers"""
        self.connect(self.ui.redWeightSlider, QtCore.SIGNAL('valueChanged(int)'), self.onRedWeightChange)
        self.connect(self.ui.greenWeightSlider, QtCore.SIGNAL('valueChanged(int)'), self.onGreenWeightChange)
        self.connect(self.ui.blueWeightSlider, QtCore.SIGNAL('valueChanged(int)'), self.onBlueWeightChange)
        
    def onRedWeightChange(self, value):
        """Red channel weight changed"""
        self._weights[0] = value / 100.
        self.ui.compositeImage.update_figure((self._weights))
        
    def onGreenWeightChange(self, value):
        """Green channel weight changed"""
        self._weights[1] = value / 100.
        self.ui.compositeImage.update_figure((self._weights))
        
    def onBlueWeightChange(self, value):
        """Blue channel weight changed"""
        self._weights[2] = value / 100.
        self.ui.compositeImage.update_figure((self._weights))     

class SunPyPlot(FigureCanvas):
    """SunPy preview image"""
    def __init__(self, map_, width, height, parent=None, dpi=100, **matplot_args):
        #self._widthHint = width
        #self._heightHint = height
        
        self._origMap = map_
        self._map = map_.resample((width, height))
        
        self.figure = self._map.plot_simple(**matplot_args)
        FigureCanvas.__init__(self, self.figure)
        
        # How can we get the canvas to preserve its aspect ratio when expanding?
        #sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        #sizePolicy.setHeightForWidth(True)
        #self.setSizePolicy(sizePolicy)
        
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        self.setSizePolicy(sizePolicy)
        self.setMinimumSize(QtCore.QSize(width, height))
        self.setMaximumSize(QtCore.QSize(width, height))   
        #FigureCanvas.updateGeometry(self)        

    def update_figure(self, weights):
        """Update RGB Composite image"""
        weightArray = np.ones((self._map.shape[0], self._map.shape[1], 3))
        
        for i in range(3):
            weightArray[:,:,i] *= weights[i]
        
        self.figure = (self._map * weightArray).astype(np.uint8).plot_simple()
        self.draw()
    
#    def heightForWidth(self, width):
#        """Preserves 1:1 aspect ratio"""
#        return width
#    
#    def sizeHint(self):
#        """Preview image default size"""
#        return QtCore.QSize(self._widthHint, self._heightHint)
    
class RGBCompositeMap(sunpy.MapCube):
    """A composite map where each color channel is associated with a separate
       datasource."""
    def __new__(cls, red, green, blue, **kwargs):
        headers = []
        data = np.zeros((red.shape[0], red.shape[1], 3), dtype=np.uint8)
        
        # convert input to maps
        for i, item in enumerate([red, green, blue]):
            if isinstance(item, BaseMap):
                map_ = item
            else:
                map_ = BaseMap.read(item)
                
            data[:,:,i] = map_
            headers.append(map_.header)

        obj = np.asarray(data).view(cls)
        obj._headers = headers

        return obj

    def __init__(self, *args, **kwargs):
        sunpy.MapCube.__init__(self, args, kwargs)
        
    def resample(self, dimensions, method='linear'):
        """Returns a new Map that has been resampled up or down
        
        See `sunpy.map.BaseMap.resample`
        """
        resampled = []
        
        for map_ in self.transpose(2, 0, 1):
            resampled.append(map_.resample(dimensions, method))

        return self.__class__(*resampled)
        
    @toggle_pylab
    def plot_simple(self, **matplot_args):
        """Plots the map object using matplotlib
        
        Parameters
        ----------
        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting the image.
        """
        fig = plt.figure(frameon=False)
        
        axes = plt.Axes(fig, [0., 0., 1., 1.])
        axes.set_axis_off()
        fig.add_axes(axes)

        axes.imshow(self, origin='lower', aspect='normal', **matplot_args)
        return fig


if __name__ == '__main__':
    sys.exit(main())