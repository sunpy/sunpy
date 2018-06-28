"""Datasource-specific classes

This is where datasource specific logic is implemented. Each mission should
have its own file with one or more classes defined.
"""
__all__ = ['CallistoSpectrogram', 'SWavesSpectrogram']

from .swaves import SWavesSpectrogram
from .callisto import CallistoSpectrogram
