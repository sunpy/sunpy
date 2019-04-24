"""
Solar Physical Models
---------------------
This module contains standard models of the Sun from various sources.

All data is saved in `astropy.table.QTable` with an added attribute:

* source : names the source of the data

Object
------
    interior : `astropy.table.QTable`
        The standard model of the solar interior.
    evolution : `astropy.table.QTable`
        The evolution as a function of time of the Sun.

References
----------
* Adapted from Turck-Chieze et al. (1988) with composition: X = 0.7046, Y = 0.2757, Z = 0.0197
"""
import astropy.units as u
from astropy.table import QTable

__all__ = ["interior", "evolution"]


# Radius -  R_sun
_radius = [0, 0.01, 0.022, 0.061, 0.090, 0.120,
           0.166, 0.202, 0.246, 0.281, 0.317,
           0.370, 0.453, 0.611, 0.7304, 0.862,
           0.965, 1.0000] * u.Rsun

# mass -  M_sun
_mass = [0, 0.0001, 0.002, 0.020, 0.057, 0.115,
         0.235, 0.341, 0.470, 0.562, 0.647, 0.748,
         0.854, 0.951, 0.9809, 0.9964, 0.9999, 1.000] * u.Msun

# luminosity -  L_sun
_luminosity = [0, 0.0009, 0.009, 0.154, 0.365,
               0.594, 0.845, 0.940, 0.985, 0.997,
               0.992, 0.9996, 1.000, 1.0, 1.0, 1.0, 1.0, 1.0] * u.Lsun

# temperature -  10^6 K
_temperature = [15.513, 15.48, 15.36, 14.404,
                13.37, 12.25, 10.53, 9.30, 8.035,
                7.214, 6.461, 5.531, 4.426, 2.981,
                2.035, 0.884, 0.1818, 0.005770] * u.MK

# density -  g cm^-3
_density = [147.74, 146.66, 142.73, 116.10, 93.35,
            72.73, 48.19, 34.28, 21.958, 15.157,
            10.157, 5.566, 2.259, 0.4483, 0.1528,
            0.042, 0.00361, 1.99e-7] * u.g*u.cm**-3

_d = {'radius': _radius, 'mass': _mass, 'luminosity': _luminosity,
      'temperature': _temperature, 'density': _density}
interior = QTable(_d)
interior.source = 'Turck-Chieze et al. (1988)'
interior.add_index('radius')

# time -  10^9 years
_time = [0, 0.143, 0.856, 1.863, 2.193, 3.020,
         3.977, 4.587, 5.506, 6.074, 6.577, 7.027,
         7.728, 8.258, 8.7566, 9.805] * u.Gyr

# luminosity -  L_sun
_tluminosity = [0.7688, 0.7248, 0.7621, 0.8156,
                0.8352, 0.8855, 0.9522, 1.0, 1.079,
                1.133, 1.186, 1.238, 1.318, 1.399,
                1.494, 1.760] * u.Lsun

# radius -  R_sun
_tradius = [0.872, 0.885, 0.902, 0.924, 0.932,
            0.953, 0.981, 1.0, 1.035, 1.059, 1.082,
            1.105, 1.143, 1.180, 1.224, 1.361] * u.Rsun

# central temperature -  10^6 K
_tcentral_temperature = [13.35, 13.46, 13.68, 14.08, 14.22, 14.60,
                         15.12, 15.51, 16.18, 16.65, 17.13, 17.62,
                         18.42, 18.74, 18.81, 19.25] * u.MK

_t = {'time': _time, 'luminosity': _tluminosity, 'radius': _tradius,
      'central temperature': _tcentral_temperature}
evolution = QTable(_t)
evolution.source = 'Unknown'
evolution.add_index('time')
