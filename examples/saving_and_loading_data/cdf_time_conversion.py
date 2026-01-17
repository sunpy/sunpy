"""
===================
Converting CDF time
===================

This example demonstrates how to convert CDF time formats to `~astropy.time.Time`.
"""
import warnings

from erfa import ErfaWarning

from sunpy.time import parse_time
from sunpy.util.exceptions import SunpyUserWarning

warnings.simplefilter('ignore', ErfaWarning)

###############################################################################
# The Common Data Format (CDF) has a number of standard time formats that are
# not supported by `astropy.time` by default. These are ``CDF_EPOCH``,
# ``CDF_EPOCH16`` and ``CDF_TT2000``.
#
# SunPy provides a way to convert these times to `~astropy.time.Time` objects
# using the `cdflib` library.

###############################################################################
# First we need to make sure that `cdflib` is installed.
try:
    import cdflib
except ImportError:
    warnings.warn("cdflib is not installed, skipping example", SunpyUserWarning)
    exit(0)

###############################################################################
# We can parse ``CDF_TT2000`` times, which are nanoseconds since J2000.
t_ns = 599572869184000000
t = parse_time(t_ns, format='cdf_tt2000')
print(f"Format: {t.format}")
print(f"Scale: {t.scale}")
print(f"Times:\n{t}")

###############################################################################
# We can also parse ``CDF_EPOCH`` times, which are milliseconds since 0000-01-01.
t_epoch = 63806836800000.0
t = parse_time(t_epoch, format='cdf_epoch')
print(f"\nCDF_EPOCH Times:\n{t}")
