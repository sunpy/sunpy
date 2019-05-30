"""
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


solar_semidiameter_angular_size = deprecated('1.0')(_sun.solar_semidiameter_angular_size)


position = deprecated('1.0')(_sun.position)


carrington_rotation_number = deprecated('1.0')(_sun.carrington_rotation_number)


true_longitude = deprecated('1.0')(_sun.true_longitude)


apparent_longitude = deprecated('1.0')(_sun.apparent_longitude)


true_latitude = deprecated('1.0')(_sun.true_latitude)


apparent_latitude = deprecated('1.0')(_sun.apparent_latitude)


mean_obliquity_of_ecliptic = deprecated('1.0')(_sun.mean_obliquity_of_ecliptic)


true_rightascension = deprecated('1.0')(_sun.true_rightascension)


true_declination = deprecated('1.0')(_sun.true_declination)


true_obliquity_of_ecliptic = deprecated('1.0')(_sun.true_obliquity_of_ecliptic)


apparent_rightascension = deprecated('1.0')(_sun.apparent_rightascension)


apparent_declination = deprecated('1.0')(_sun.apparent_declination)


print_params = deprecated('1.0')(_sun.print_params)
