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


solar_semidiameter_angular_size = deprecated('1.0', name='solar_semidiameter_angular_size',
                                             alternative='sunpy.coordinates.sun.angular_radius'
                                            )(_sun.angular_radius)


position = deprecated('1.0', name='position',
                      alternative='sunpy.coordinates.sun.sky_position'
                     )(_sun.sky_position)


carrington_rotation_number = deprecated('1.0',
                                        alternative='sunpy.coordinates.sun.carrington_rotation_number'
                                       )(_sun.carrington_rotation_number)


true_longitude = deprecated('1.0',
                            alternative='sunpy.coordinates.sun.true_longitude'
                           )(_sun.true_longitude)


apparent_longitude = deprecated('1.0',
                                alternative='sunpy.coordinates.sun.apparent_longitude'
                               )(_sun.apparent_longitude)


true_latitude = deprecated('1.0',
                           alternative='sunpy.coordinates.sun.true_latitude'
                          )(_sun.true_latitude)


apparent_latitude = deprecated('1.0',
                               alternative='sunpy.coordinates.sun.apparent_latitude'
                              )(_sun.apparent_latitude)


apparent_latitude = deprecated('1.0',
                               alternative='sunpy.coordinates.sun.apparent_latitude'
                              )(_sun.apparent_latitude)


mean_obliquity_of_ecliptic = deprecated('1.0',
                                        alternative='sunpy.coordinates.sun.mean_obliquity_of_ecliptic'
                                       )(_sun.mean_obliquity_of_ecliptic)


true_rightascension = deprecated('1.0',
                                 alternative='sunpy.coordinates.sun.true_rightascension'
                                )(_sun.true_rightascension)


true_declination = deprecated('1.0',
                              alternative='sunpy.coordinates.sun.true_declination'
                             )(_sun.true_declination)


true_obliquity_of_ecliptic = deprecated('1.0',
                                        alternative='sunpy.coordinates.sun.true_obliquity_of_ecliptic'
                                       )(_sun.true_obliquity_of_ecliptic)


apparent_rightascension = deprecated('1.0',
                                     alternative='sunpy.coordinates.sun.apparent_rightascension'
                                    )(_sun.apparent_rightascension)


apparent_declination = deprecated('1.0',
                                  alternative='sunpy.coordinates.sun.apparent_declination'
                                 )(_sun.apparent_declination)


print_params = deprecated('1.0',
                          alternative='sunpy.coordinates.sun.print_params'
                         )(_sun.print_params)
