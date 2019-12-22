from .dataretriever.attrs import goes
from .jsoc import attrs as jsoc
from .vso import attrs as vso
from .vso.attrs import Detector, Instrument, Level, Resolution, Sample, Time, Wavelength

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'Sample', 'Detector', 'Resolution', 'vso',
           'jsoc', 'goes']
