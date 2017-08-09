"""
The Interface Region Imaging Spectrometer (IRIS) instrument routines.

Currently this module only includes a simple hack to convert SJI mode files to
a `sunpy.map.MapCube` instance.

.. note::

    More comprehensive IRIS tools are now being developed in the `IRISPy
    <https://github.com/sunpy/irispy>`__ affiliated pacakge.
"""
from .iris import *
