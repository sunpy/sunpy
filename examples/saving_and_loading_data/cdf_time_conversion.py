
"""
============================================
Converting CDF time formats to Astropy Times
============================================

This example demonstrates how to convert CDF time formats (CDF_EPOCH, CDF_EPOCH16, TT2000)
to `astropy.time.Time` objects using `~sunpy.time.parse_time`.
This is useful when working with CDF files that store time in these specific formats.
"""
import numpy as np

import sunpy.time

################################################################################
# The Common Data Format (CDF) is widely used in space physics.
# It defines specific time formats like CDF_EPOCH (milliseconds since 0 AD),
# CDF_EPOCH16 (picoseconds since 0 AD), and TT2000 (nanoseconds since J2000).
#
# Often these are read as arrays of numbers. SunPy can convert these into
# `astropy.time.Time` objects if ``cdflib`` is installed.

################################################################################
# Example with CDF_TT2000 (nanoseconds since J2000)
# Let's create a numpy array representing TT2000 times.
# 0 corresponds to J2000 (2000-01-01 12:00:00 TT).
tt2000_times = np.array([0, 1000000000, 60000000000], dtype=np.int64)

# We can parse these using `parse_time` by specifying the format 'cdf_tt2000'.
# Note: This requires the optional dependency 'cdflib'.
try:
    time_obj = sunpy.time.parse_time(tt2000_times, format='cdf_tt2000')

    print(f"Format: {time_obj.format}")
    print(f"Scale: {time_obj.scale}")
    print("Times:")
    print(time_obj.iso)
except ImportError as e:
    print(f"Could not parse CDF times: {e}")

################################################################################
# Example with CDF_EPOCH (milliseconds since 0 AD)
# 6.383e13 milliseconds is roughly year 2023.
cdf_epoch_times = np.array([63830592000000.0])

try:
    time_epoch = sunpy.time.parse_time(cdf_epoch_times, format='cdf_epoch')
    print("\nCDF_EPOCH Times:")
    print(time_epoch.iso)
except ImportError:
    pass

################################################################################
# These `astropy.time.Time` objects can now be used for analysis or plotting
# with other SunPy and Astropy tools.
