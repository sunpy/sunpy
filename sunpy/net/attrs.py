# -*- coding: utf-8 -*-
"""
This module is a wrapper to all the query attributes that can be used with the `~sunpy.net.dataretriever.Fido` interface.

There are four attributes that are global to all the clients:

Global Attributes
^^^^^^^^^^^^^^^^^

.. autosummary::

    Time
    Instrument
    Wavelength
    Level

"""

from .vso import attrs as vso
from .jsoc import attrs as jsoc

from .vso.attrs import Time, Instrument, Wavelength, Level

__all__ = ['Time', 'Instrument', 'Wavelength', 'Level', 'vso', 'jsoc']
