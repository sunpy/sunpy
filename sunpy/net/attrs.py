# -*- coding: utf-8 -*-

from .vso import attrs as vso
from .jsoc import attrs as jsoc

from .vso.attrs import Time, Instrument, Wavelength, Level, Source, Physobs, Detector, Filter

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'vso', 'jsoc',
           'Source', 'Physobs', 'Detector', 'Filter']
