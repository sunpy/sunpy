#!/usr/bin/env python
#-*- coding:utf-8 -*-
#pylint: disable=W0223
"""
RGB Composite Image Demo

Last Update: Feb 13, 2012
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
    * Wavelength adjustment support.
    * Refactor/simplify
    * Make all images expand to fill available space.
"""
import sys
import datetime
import sunpy
import Ui_RGBComposite
import numpy as np
from sunpy.net.helioviewer import HelioviewerClient
import sunpy.map
from sunpy.map import Map
from sunpy.util import toggle_pylab
from PyQt4 import QtGui, QtCore
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

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

        # Helioviewer client
        self._hv = HelioviewerClient()
        self._datasources = None

        # Loaded images
        self.red = None
        self.green = None
        self.blue = None

        # Color channel weights
        self._weights=[1., 1., 1.]

        # Setup UI
        self._load_data_sources()

        # Load initial data
        self._load_defaults()

        # Initialize event handlers
        self._initEvents()

    def _load_data_sources(self):
        """Downloads and displays latest images for default wavelengths"""
        self._datasources = self._hv.get_data_sources()['SDO']['AIA']['AIA']
        sorted_datasources = sorted(self._datasources, key=int)

        for wl in sorted_datasources:
            self.ui.redWavelengthSelect.addItem(wl, self._datasources[wl])
            self.ui.greenWavelengthSelect.addItem(wl, self._datasources[wl])
            self.ui.blueWavelengthSelect.addItem(wl, self._datasources[wl])

        # Default wavelengths: 304, 193, 171
        self.ui.redWavelengthSelect.setCurrentIndex(5)
        self.ui.greenWavelengthSelect.setCurrentIndex(3)
        self.ui.blueWavelengthSelect.setCurrentIndex(2)

    def _load_defaults(self):
        """Load initial images"""
        now = datetime.datetime.utcnow()
        self.ui.dateTimeEdit.setDateTime(now)

        r = self._hv.download_jp2(now, sourceId=self._datasources['304']['sourceId'])
        g = self._hv.download_jp2(now, sourceId=self._datasources['193']['sourceId'])
        b = self._hv.download_jp2(now, sourceId=self._datasources['171']['sourceId'])

        self.red = sunpy.map.Map(r)
        self.green = sunpy.map.Map(g)
        self.blue = sunpy.map.Map(b)

        self._updateRedPreview()
        self._updateGreenPreview()
        self._updateBluePreview()

        self._createCompositeImage()

    def _updateCompositeImage(self):
        """Updates the composite image"""
        # Remove old composite
        rgb = RGBCompositeMap(self.red, self.green, self.blue)
        self.ui.compositeContainer.removeWidget(self.ui.compositeImage)
        self.ui.compositeImage.close()

        # Plot new one
        self.ui.compositeImage = RGBCompositePlot(rgb, 512, 512)
        self.ui.compositeImage.set_rgb_weights(self._weights)
        self.ui.compositeContainer.addWidget(self.ui.compositeImage, 1)
        self.ui.compositeContainer.update()

    def _createCompositeImage(self):
        """Creates an initial composite image plot"""
        rgb = RGBCompositeMap(self.red, self.green, self.blue)
        self.ui.compositeContainer.removeWidget(self.ui.compositePlaceholder)
        self.ui.compositePlaceholder.close()
        self.ui.compositeImage = RGBCompositePlot(rgb, 512, 512)
        self.ui.compositeContainer.addWidget(self.ui.compositeImage, 1)
        self.ui.compositeContainer.update()

    def _updateRedPreview(self):
        """Updates the red preview image"""
        if hasattr(self.ui, "redPreviewImage"):
            self.ui.redPreview.removeWidget(self.ui.redPreviewImage)
            self.ui.redPreviewImage.close()
        else:
            self.ui.redPreview.removeWidget(self.ui.redPlaceholder)
            self.ui.redPlaceholder.close()

        self.ui.redPreviewImage = SunPyPlot(self.red, 256, 256) #cmap=cm.Reds_r
        self.ui.redPreview.addWidget(self.ui.redPreviewImage, 1)
        self.ui.redPreview.update()

    def _updateGreenPreview(self):
        """Updates the green preview image"""
        if hasattr(self.ui, "greenPreviewImage"):
            self.ui.greenPreview.removeWidget(self.ui.greenPreviewImage)
            self.ui.greenPreviewImage.close()
        else:
            self.ui.greenPreview.removeWidget(self.ui.greenPlaceholder)
            self.ui.greenPlaceholder.close()

        self.ui.greenPreviewImage = SunPyPlot(self.green, 256, 256) #cmap=cm.Greens_r
        self.ui.greenPreview.addWidget(self.ui.greenPreviewImage, 1)
        self.ui.greenPreview.update()

    def _updateBluePreview(self):
        """Updates the blue preview image"""
        if hasattr(self.ui, "bluePreviewImage"):
            self.ui.bluePreview.removeWidget(self.ui.bluePreviewImage)
            self.ui.bluePreviewImage.close()
        else:
            self.ui.bluePreview.removeWidget(self.ui.bluePlaceholder)
            self.ui.bluePlaceholder.close()

        self.ui.bluePreviewImage = SunPyPlot(self.blue, 256, 256) #cmap=cm.Blues_r
        self.ui.bluePreview.addWidget(self.ui.bluePreviewImage, 1)
        self.ui.bluePreview.update()

    def _initEvents(self):
        """Initialize event handlers"""
        self.connect(self.ui.redWeightSlider, QtCore.SIGNAL('valueChanged(int)'), self.onRedWeightChange)
        self.connect(self.ui.greenWeightSlider, QtCore.SIGNAL('valueChanged(int)'), self.onGreenWeightChange)
        self.connect(self.ui.blueWeightSlider, QtCore.SIGNAL('valueChanged(int)'), self.onBlueWeightChange)
        self.connect(self.ui.dateTimeEdit, QtCore.SIGNAL('dateTimeChanged(QDateTime)'), self.onDateChange)
        self.connect(self.ui.dateTimeEdit, QtCore.SIGNAL('dateTimeChanged(QDateTime)'), self.onDateChange)
        self.connect(self.ui.actionSave, QtCore.SIGNAL('activated()'), self.onSaveClick)

    def onRedWeightChange(self, value):
        """Red channel weight changed"""
        self._weights[0] = value / 100.
        self.ui.compositeImage.set_rgb_weights((self._weights))

    def onGreenWeightChange(self, value):
        """Green channel weight changed"""
        self._weights[1] = value / 100.
        self.ui.compositeImage.set_rgb_weights((self._weights))

    def onBlueWeightChange(self, value):
        """Blue channel weight changed"""
        self._weights[2] = value / 100.
        self.ui.compositeImage.set_rgb_weights((self._weights))

    def onDateChange(self, qdatetime):
        """Updates the images when the date is changed"""
        dt = qdatetime.toPyDateTime()

        r = self._hv.download_jp2(dt, sourceId=self._datasources['304']['sourceId'])
        g = self._hv.download_jp2(dt, sourceId=self._datasources['193']['sourceId'])
        b = self._hv.download_jp2(dt, sourceId=self._datasources['171']['sourceId'])

        self.red = sunpy.map.Map(r)
        self.green = sunpy.map.Map(g)
        self.blue = sunpy.map.Map(b)

        self._updateRedPreview()
        self._updateGreenPreview()
        self._updateBluePreview()

        self._updateCompositeImage()

    def onSaveClick(self):
        """Save the composite image"""
        filename = QtGui.QFileDialog.getSaveFileName(self, "Save image", "composite.png")
        self.ui.compositeImage.figure.savefig(str(filename))

class SunPyPlot(FigureCanvas):
    """SunPy preview image"""
    def __init__(self, map_, width, height, parent=None, dpi=100,
                 **matplot_args): #pylint: disable=W0613
        #self._widthHint = width
        #self._heightHint = height

        self._origMap = map_
        self._map = map_.resample((width, height))

        # Old way (segfaults in some environements)
        #self.figure = self._map.plot_simple(**matplot_args)
        #FigureCanvas.__init__(self, self.figure)

        self.figure = Figure()
        self._map.plot(figure=self.figure, basic_plot=True, **matplot_args)
        self.axes = self.figure.gca()
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

class PreviewImagePlot(SunPyPlot):
    """Qt representation of a preview thumbnail"""
    def __init__(self, map_, width, height, parent=None, dpi=100, **matplot_args):
        SunPyPlot.__init__(self, map_, width, height, parent=None, dpi=100, **matplot_args)

class RGBCompositePlot(SunPyPlot):
    """Qt representation of an RGB composite image"""
    def __init__(self, map_, width, height, parent=None, dpi=100, **matplot_args):
        SunPyPlot.__init__(self, map_, width, height, parent=None, dpi=100, **matplot_args)

    def set_rgb_weights(self, weights):
        """Update RGB Composite image weights"""
        weightArray = np.ones((self._map.shape[0], self._map.shape[1], 3))

        # Rescale data using new weights
        for i in range(3):
            weightArray[:,:,i] *= weights[i]
        new_data = (self._map * weightArray).astype(np.uint8)

        # Update AxesImage data (faster than calling imshow again)
        ax = self.figure.get_axes()[0].images[0]
        ax.set_data(new_data)
        self.draw()

#    def heightForWidth(self, width):
#        """Preserves 1:1 aspect ratio"""
#        return width
#
#    def sizeHint(self):
#        """Preview image default size"""
#        return QtCore.QSize(self._widthHint, self._heightHint)

class RGBCompositeMap(sunpy.map.MapCube):
    """A composite map where each color channel is associated with a separate
       datasource."""
    def __new__(cls, red, green, blue, **kwargs):
        headers = []
        data = np.zeros((red.shape[0], red.shape[1], 3), dtype=np.uint8)

        # convert input to maps
        for i, item in enumerate([red, green, blue]):
            if isinstance(item, Map):
                map_ = item
            else:
                map_ = Map.read(item)

            data[:,:,i] = map_
            headers.append(map_.get_header(original=True))

        obj = np.asarray(data).view(cls)
        obj._headers = headers

        return obj

    def __init__(self, *args, **kwargs):
        sunpy.MapCube.__init__(self, args, kwargs)

    def resample(self, dimensions, method='linear'):
        """Returns a new Map that has been resampled up or down

        See `sunpy.map.Map.resample`
        """
        resampled = []

        for map_ in self.transpose(2, 0, 1): #pylint: disable=E1101
            resampled.append(map_.resample(dimensions, method))

        return self.__class__(*resampled)

    @toggle_pylab
    def plot(self, figure=None, basic_plot=None, **matplot_args):
        """Plots the map object using matplotlib

        Parameters
        ----------
        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting.
        """
        if figure is None:
            figure = plt.figure(frameon=False)

        axes = plt.Axes(figure, [0., 0., 1., 1.])
        axes.set_axis_off()
        figure.add_axes(axes)

        axes.imshow(self, origin='lower', aspect='normal', **matplot_args)
        return figure

if __name__ == '__main__':
    sys.exit(main())
