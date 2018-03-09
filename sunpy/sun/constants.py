"""
Fundamental Solar Physical Constants
------------------------------------
These constants are taken from various sources. The structure of this module is heavily
based on if not directly copied from the SciPy constants module but contains Solar
Physical constants.
"""

from __future__ import absolute_import, division, print_function
from astropy.table import Table
from sunpy.extern.six import iteritems

from sunpy.sun import _constants as _con  # pylint: disable=E0611

__all__ = [
    'get', 'find', 'print_all', 'spectral_classification', 'au', 'mass', 'equatorial_radius',
    'volume', 'surface_area', 'average_density', 'equatorial_surface_gravity',
    'effective_temperature', 'luminosity', 'mass_conversion_rate', 'escape_velocity', 'sfu',
    'average_angular_size'
]

constants = _con.physical_constants


def get(key):
    """
    Retrieve a constant by key. This is just a short cut into a dictionary.

    Parameters
    ----------
    key : Python string or unicode
        Key in dictionary in `constants`

    Returns
    -------
    constant :  `~astropy.units.Constant`

    See Also
    --------
    _constants : Contains the description of `constants`, which, as a
        dictionary literal object, does not itself possess a docstring.

    Examples
    --------
    >>> from sunpy.sun import constants
    >>> constants.get('mass')
    <<class 'astropy.constants.iau2012.IAU2012'> name='Solar mass' value=1.9891e+30 uncertainty=5e+25 unit='kg' reference="Allen's Astrophysical Quantities 4th Ed.">
    """
    return constants[key]


def find(sub=None):
    """
    Return list of constants keys containing a given string

    Parameters
    ----------
    sub : str, unicode
        Sub-string to search keys for.  By default, return all keys.

    Returns
    -------
    keys : None or list

    See Also
    --------
    _constants : Contains the description of `constants`, which, as a
        dictionary literal object, does not itself possess a docstring.

    """
    if sub is None:
        result = list(constants.keys())
    else:
        result = [key for key in constants \
                 if sub.lower() in key.lower()]

    result.sort()
    return result


def print_all(key=None):
    """
    Provides a table of the complete list of constants.

    Parameters
    ----------
    key : Python string or unicode
        Key in dictionary `constants`

    Returns
    -------
    table : `astropy.table.Table`
    """
    data_rows = []
    for key, this_constant in iteritems(constants):
        data_rows.append([
            key, this_constant.name, this_constant.value, this_constant.uncertainty,
            str(this_constant.unit), this_constant.reference
        ])

    t = Table(rows=data_rows, names=('key', 'name', 'value', 'uncertainty', 'unit', 'Reference'))
    return t


# Spectral class is not included in physical constants since it is not a number
spectral_classification = 'G2V'

au = astronomical_unit = get('mean distance')

# The following variables from _gets are brought out by making them
# accessible through a call such as sun.volume
mass = get('mass')
equatorial_radius = radius = get('radius')
volume = get('volume')
surface_area = get('surface area')
average_density = density = get('average density')
equatorial_surface_gravity = surface_gravity = get('surface gravity')
effective_temperature = get('effective temperature')
luminosity = get('luminosity')
mass_conversion_rate = get('mass conversion rate')
escape_velocity = get('escape velocity')

sfu = get('solar flux unit')

# Observable parameters
average_angular_size = get('average angular size')
