"""
Color options widget

Provides color map selector, clipping and scaling options,
and histogram.

Author: Matt Earnshaw <matt@earnshaw.org.uk>
"""

import sunpy.cm
import matplotlib.cm
import ui_coloroptions
import matplotlib.colors
from PyQt4.QtCore import pyqtSignature
from PyQt4.QtGui import QWidget
from matplotlib import pyplot as plt


class ColorOptions(QWidget, ui_coloroptions.Ui_ColorOptions):
    
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setupUi(self)
        
        self.params = {}
        self.vmin = self.clipMinDoubleSpinBox.value()
        self.vmax = self.clipMaxDoubleSpinBox.value()
        self.scaling_type = self.scalingComboBox.currentText()

        # Populate list widget with SunPy colormaps
        for cmap in sunpy.cm.cmlist:
            self.cmListWidget.addItem(cmap)

        # Populate list widget with MPL colormaps
        self.mpl_cmaps = sorted(m for m in plt.cm.datad if not m.endswith("_r"))
        for cmap in self.mpl_cmaps:
            self.cmListWidget.addItem(cmap)


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

        self.update_plot(cmap=new_cmap)

    @pyqtSignature("double")
    def on_clipMinDoubleSpinBox_valueChanged(self, vmin):
        self.vmax = self.clipMaxDoubleSpinBox.value()
        self.vmin = self.clipMinDoubleSpinBox.value()
        self.normalize_plot()

    @pyqtSignature("double")
    def on_clipMaxDoubleSpinBox_valueChanged(self, vmax):
        self.vmin = self.clipMinDoubleSpinBox.value()
        self.vmax = self.clipMaxDoubleSpinBox.value()
        self.normalize_plot()

    @pyqtSignature("QString")
    def on_scalingComboBox_currentIndexChanged(self, scaling_type):
        self.scaling_type = str(scaling_type)
        self.normalize_plot()
        
    def normalize_plot(self):
        # This is not robust! Ideas?
        if self.vmin < self.vmax:
            if self.scaling_type == "Linear":
                norm = matplotlib.colors.Normalize(vmin=self.vmin, vmax=self.vmax)
                self.update_plot(norm=norm)
            elif self.scaling_type == "Logarithmic":
                norm = matplotlib.colors.LogNorm(vmin=self.vmin, vmax=self.vmax)
                self.update_plot(norm=norm)
        else:
            self.clipMinDoubleSpinBox.setValue(self.vmax)

    def update_plot(self, **kwargs):
        """ Provide interface to plot method of current tab's map """
        current_tab = self.current_tab()
        current_tab.canvas.figure = current_tab._map.plot(**kwargs)
        current_tab.canvas.draw()

    def current_tab(self):
        return self.window().tabWidget.currentWidget()
