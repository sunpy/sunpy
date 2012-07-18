# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import numpy as np

class Spectrum(np.ndarray):
    def __new__(cls, data, freq_axis):
        return np.asarray(data).view(cls)

    def __init__(self, data, freq_axis):
    	self.freq_axis = freq_axis
