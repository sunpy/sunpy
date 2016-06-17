# -*- coding: utf-8 -*-

from .vso import attrs as vso
from .jsoc import attrs as jsoc

from .vso.attrs import Time, Instrument, Wavelength, Level, Physobs, Source, Detector, Filter

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'vso', 'jsoc',
           'Physobs', 'Source', 'Detector', 'Filter']
