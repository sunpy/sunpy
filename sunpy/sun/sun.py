"""
.. deprecated:: 1.0
    This module is deprecated and may be removed in a future version.
        Use sunpy.coordinates.sun instead.

This module provides Sun-related parameters.
"""
from sunpy.coordinates import sun as _sun
from sunpy.util.decorators import deprecated

__all__ = [
    "solar_semidiameter_angular_size", "position", "carrington_rotation_number",
    "true_longitude", "apparent_longitude", "true_latitude", "apparent_latitude",
    "mean_obliquity_of_ecliptic", "true_rightascension", "true_declination",
    "true_obliquity_of_ecliptic", "apparent_rightascension", "apparent_declination",
    "print_params"
]


# The names for the functions in sunpy.coordinates.sun
_new_module = 'sunpy.coordinates.sun.'
_new_names = [
    "angular_radius", "sky_position", "carrington_rotation_number",
    "true_longitude", "apparent_longitude", "true_latitude", "apparent_latitude",
    "mean_obliquity_of_ecliptic", "true_rightascension", "true_declination",
    "true_obliquity_of_ecliptic", "apparent_rightascension", "apparent_declination",
    "print_params"
]

# Create a deprecation hook for each of the functions
for old, new in zip(__all__, _new_names):
    vars()[old] = deprecated('1.0', name=old, alternative=_new_module + new)(getattr(_sun, new))
    vars()[old].__module__ = __name__  # so that docs think that the function is local
