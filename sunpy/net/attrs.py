# -*- coding: utf-8 -*-

from .vso import attrs as vso
from .jsoc import attrs as jsoc

from .vso.attrs import Time, Instrument, Wavelength

__all__ = ['Time', 'Instrument', 'Wavelength', 'vso', 'jsoc']
