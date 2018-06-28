# -*- coding: utf-8 -*-

from .vso import attrs as vso
from .jsoc import attrs as jsoc
from .vso.attrs import Time, Level, Sample, Instrument, Wavelength
from .dataretriever.attrs import goes

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'Sample', 'vso', 'jsoc', 'goes']
