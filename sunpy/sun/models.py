"""
Solar Physical Models
---------------------
This module contains models of the Sun from various sources:

* ``interior``: `~astropy.table.QTable` of the structure of the solar interior
  as defined in Table 7 of :cite:t:`turck-chieze_revisiting_1988`.
* ``evolution``: `~astropy.table.QTable` of the evolution of the Sun over time
  as defined in Table 6 of :cite:t:`turck-chieze_revisiting_1988`.
* :func:`~sunpy.sun.models.differential_rotation`: Function for calculating
  solar differential rotation for different models
"""
import numpy as np

import astropy.units as u
from astropy.coordinates import Longitude
from astropy.table import QTable

from sunpy.sun.constants import sidereal_rotation_rate

__all__ = ["interior", "evolution", "differential_rotation"]


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
evolution.source = 'Turck-Chieze et al. (1988)'
evolution.add_index('time')


@u.quantity_input
def differential_rotation(duration: u.s, latitude: u.deg, *, model='howard', frame_time='sidereal'):
    r"""
    Computes the change in longitude over a duration for a given latitude.

    Since the Sun is not a rigid body, different heliographic latitudes rotate with
    different periods. This is known as solar differential rotation.

    Parameters
    ----------
    duration : `~astropy.units.Quantity`
        Amount of time to rotate over.
    latitude : `~astropy.units.Quantity`
        Heliographic latitude.
    model : `str`
        The differential-rotation model to use.

        One of:

        | ``howard`` : Use values from :cite:t:`howard_solar_1990`
        | ``snodgrass`` : Use values from :cite:t:`snodgrass_magnetic_1983`
        | ``allen`` : Use values from Allen's Astrophysical Quantities, and simpler equation.
        | ``rigid`` : Use values from `~sunpy.sun.constants.sidereal_rotation_rate`.

    frame_time : `str`
        If ``'sidereal'``, returns the change in longitude as referenced to distant
        stars. If ``'synodic'``, returns the apparent change in longitude as
        observed by the average orbital motion of Earth, which results in a slower
        rotation rate. Defaults to ``'sidereal'``.

    Returns
    -------
    `~astropy.units.Quantity`
        The change in longitude

    Notes
    -----
    The rotation rate at a heliographic latitude :math:`\theta` is given by

    .. math::

        A + B \sin^{2} \left (\theta \right ) + C \sin^{4} \left ( \theta \right )

    where :math:`A, B, C` are constants that depend on the model:

    ========= ======= ====== ====== ==========
    Model     A       B      C      Unit
    ========= ======= ====== ====== ==========
    howard    2.894   -0.428 -0.370 microrad/s
    snodgrass 2.851   -0.343 -0.474 microrad/s
    allen     14.44   -3.0   0      deg/day
    rigid     14.1844 0      0      deg/day
    ========= ======= ====== ====== ==========

    1 microrad/s is approximately 4.95 deg/day.
    See also the comparisons in :cite:t:`beck_comparison_2000`.

    Examples
    --------
    .. minigallery:: sunpy.sun.models.differential_rotation

    Default rotation calculation over two days at 30 degrees latitude:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from sunpy.sun.models import differential_rotation
    >>> differential_rotation(2 * u.day, 30 * u.deg)
    <Longitude 27.36432679 deg>

    Default rotation over two days for a number of latitudes:

    >>> differential_rotation(2 * u.day, np.linspace(-70, 70, 20) * u.deg)
    <Longitude [22.05449682, 23.03214991, 24.12033958, 25.210281  ,
                26.21032832, 27.05716463, 27.71932645, 28.19299667,
                28.49196765, 28.63509765, 28.63509765, 28.49196765,
                28.19299667, 27.71932645, 27.05716463, 26.21032832,
                25.210281  , 24.12033958, 23.03214991, 22.05449682] deg>

    With rotation model 'allen':

    >>> differential_rotation(2 * u.day, np.linspace(-70, 70, 20) * u.deg, model='allen')
    <Longitude [23.58186667, 24.14800185, 24.82808733, 25.57737945,
            26.34658134, 27.08508627, 27.74430709, 28.28087284,
            28.6594822 , 28.85522599, 28.85522599, 28.6594822 ,
            28.28087284, 27.74430709, 27.08508627, 26.34658134,
            25.57737945, 24.82808733, 24.14800185, 23.58186667] deg>
    """

    latitude = latitude.to(u.deg)

    sin2l = (np.sin(latitude))**2
    sin4l = sin2l**2

    rot_params = {'howard': [2.894, -0.428, -0.370] * u.urad / u.second,
                  'snodgrass': [2.851, -0.343, -0.474] * u.urad / u.second,
                  'allen': [14.44, -3.0, 0] * u.deg / u.day,
                  'rigid': sidereal_rotation_rate * [1, 0, 0]
                  }

    if model not in rot_params:
        raise ValueError("model must equal one of "
                         "{{ {} }}".format(" , ".join(rot_params.keys())))

    A, B, C = rot_params[model]

    # This calculation of the rotation assumes a sidereal frame time.
    rotation = (A + B * sin2l + C * sin4l) * duration

    # Applying this correction assumes that the observer is on the Earth,
    # and that the Earth is at the same distance from the Sun at all times
    # during the year.
    if frame_time == 'synodic':
        rotation -= 0.9856 * u.deg / u.day * duration

    return Longitude(rotation.to(u.deg))
