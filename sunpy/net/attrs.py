# -*- coding: utf-8 -*-

from .vso import attrs as vso
from .jsoc import attrs as jsoc
from .dataretriever.attrs import goes

from .vso.attrs import Time, Instrument, Wavelength, Level, Sample, Detector, Resolution
from .vso.attrs import Source, Physobs

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'Sample', 'Detector', 'Resolution', 'vso',
           'jsoc', 'goes', 'Source', 'Physobs']
