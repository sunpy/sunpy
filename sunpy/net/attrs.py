# -*- coding: utf-8 -*-

from .vso import attrs as vso
from .jsoc import attrs as jsoc

from .vso.attrs import Time, Instrument, Wavelength, Level, Source, Detector, Physobs

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'vso', 'jsoc',
           'Source', 'Detector', 'Physobs']
