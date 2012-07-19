# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import numpy as np

from matplotlib import pyplot as plt

class Spectrum(np.ndarray):
    def __new__(cls, data, freq_axis):
        return np.asarray(data).view(cls)

    def __init__(self, data, freq_axis):
    	self.freq_axis = freq_axis

    def plot(self, overlays=[], colorbar=True, **matplotlib_args):
        # [] as default argument is okay here because it is only read.
        # pylint: disable=W0102,R0914

        figure = plt.figure(frameon=True)
        axes = figure.add_subplot(111)
        
        params = {}

        params.update(matplotlib_args)

        axes.plot(self.freq_axis, self, **params)

        for overlay in overlays:
            figure, axes = overlay(figure, axes)
        return figure
