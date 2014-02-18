"""
Contains astronomical and physical constants for use in SunPy or other
places.

The package contains a `~sunpy.sun._cgs` and `~sunpy.sun._si`
module that define constants in CGS and SI units, respectively.  
A typical use case might be::

    from sunpy.sun._cgs import physical_constants as con

"""

# Not sure why this is necessary right now, but SunPy doesn't find
# the CGS module otherwise.
# TODO: Note that with the new _cgs and _si packages, the use of constants.py
# which was copied over from SciPy does not support using the cgs package as
# it is hard coded to import the si units. Should we move away from this 
# paradigm?

from __future__ import absolute_import

from sunpy.sun.sun import *

from . import _si
from . import constants

