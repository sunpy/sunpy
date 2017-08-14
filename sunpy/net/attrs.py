# -*- coding: utf-8 -*-

from .vso import attrs as vso
from .jsoc import attrs as jsoc
from .dataretriever.attrs import goes

from .vso.attrs import Time, Instrument, Wavelength, Level

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'vso', 'jsoc', 'goes']
