# -*- coding: utf-8 -*-

from .vso import attrs as vso
from .jsoc import attrs as jsoc
from .dataretriever.attrs import goes, gbm

from .vso.attrs import Time, Instrument, Wavelength, Level, Sample, Detector

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'Sample', 'Detector', 'vso', 'jsoc',
           'goes', 'gbm']
