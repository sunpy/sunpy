# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
"""
A spectral cube is, at its most basic, a 2D array of Spectrum objects with an
associated coordinate system. The name cube is a bit misleading beacuse the
shape of the structure isn't necessarily a cuboid - each spectrum may have
different wavelength axes.
"""

import gwcs
from sunpy.spectra.spectrum import Spectrum


class SpectralCube():
    def __init__(self, spectra, wcs):
        self.spectra = spectra
        self.wcs = wcs
