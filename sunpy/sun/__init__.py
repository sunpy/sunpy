"""
Contains astronomical and physical constants for use in SunPy or other
places.

The package contains a `~sunpy.sun._constants`
module that define the constants.
A typical use case might be::

    from sunpy.sun._constants import physical_constants as con

"""

from __future__ import absolute_import

from sunpy.sun.sun import *

from . import _constants
from . import constants
