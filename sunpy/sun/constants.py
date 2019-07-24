"""
This module provides fundamental solar physical constants.
"""
from astropy.table import Table

from sunpy.sun import _constants as _con

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
    key : `str`
        Key in dictionary in ``constants``.

    Returns
    -------
    constant : `~astropy.units.Constant`

    See Also
    --------
    `sunpy.sun.constants`
        Contains the description of ``constants``, which, as a dictionary literal object, does not itself possess a docstring.

    Examples
    --------
    >>> from sunpy.sun import constants
    >>> constants.get('mass')
    <<class 'astropy.constants.iau2015.IAU2015'> name='Solar mass' value=1.9884754153381438e+30 uncertainty=9.236140093538353e+25 unit='kg' reference='IAU 2015 Resolution B 3 + CODATA 2014'>
    """
    return constants[key]


def find(sub=None):
    """
    Return list of constants keys containing a given string.

    Parameters
    ----------
    sub : `str`, optional
        Sub-string to search keys for. By default set to `None` and returns all keys.

    Returns
    -------
    `None`, `list`
        The matching keys.

    See Also
    --------
    `sunpy.sun.constants`
        Contains the description of ``constants``, which, as a dictionary literal object, does not itself possess a docstring.
    """
    if sub is None:
        result = list(constants.keys())
    else:
        result = [key for key in constants if sub.lower() in key.lower()]

    result.sort()
    return result


def print_all():
    """
    Provides a table of the complete list of constants.

    Returns
    -------
    `astropy.table.Table`
    """
    data_rows = []
    for key, this_constant in constants.items():
        data_rows.append([
            key, this_constant.name, this_constant.value, this_constant.uncertainty,
            str(this_constant.unit), this_constant.reference
        ])

    t = Table(rows=data_rows, names=('key', 'name', 'value', 'uncertainty', 'unit', 'Reference'))
    return t


# Add a list of constants to the docs
_lines = [
    'The following constants are available:\n',
    '====================== ============== ================ =================================',
    '         Name              Value            Unit                 Description',
    '====================== ============== ================ =================================',
]
for key, const in constants.items():
    _lines.append('{0:^22} {1:^14.9g} {2:^16} {3}'.format(
        key, const.value, const._unit_string, const.name))
_lines.append(_lines[1])
if __doc__ is not None:
    __doc__ += '\n'.join(_lines)

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
