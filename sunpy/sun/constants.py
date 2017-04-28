"""
Fundamental Solar Physical Constants
------------------------------------
These constants are taken from various sources. The structure of this module is heavily
based on if not directly copied from the SciPy constants module but contains Solar
Physical constants.

Object
------
    constants : dict
        A dictionary containing physical constants. Keys are the names
        of physical constants, values are tuples (value, units, uncertainty). The dictionary
        contains the following solar physical constants:

    average density:
         The average density of the Sun.
    average_angular_size:
        The average angular size of the Sun as seen from Earth in arcseconds.
    effective temperature:
        The effective black-body temperature of the Sun in Kelvin.
    oblateness:
        The ellipticity of the Sun.
    escape velocity:
        The velocity which an object needs to escape from the gravitational pull of the Sun.
    luminosity:
        The luminosity of the Sun.
    mass:
        The mass of the Sun.
    mass conversion rate:
        The rate at which the Sun converts mass to energy.
    mean energy production:
        The mean rate at which the Sun produces energy.
    mean intensity:
        The mean intensity of the Sun.
    metallicity:
        The metallicity of the Sun.
    radius:
        The radius of the Sun at the equator.
    solar flux unit:
        The definition of a solar flux unit.
    sunspot cycle:
        The average duration of the solar activity cycle.
    surface area:
        The surface area of the Sun.
    surface gravity:
        The gravitational acceleration at the surface of the Sun as measured at the equator.
    visual magnitude:
       A measure of the Sun's brightness as seen by an observer on Earth without the
       presence of the atmosphere.
    volume:
        The volume of the Sun.

Attributes
----------
A number of variables from constants are made available for convenience as
attributes.

Websites
--------
| http://books.google.com/books?id=4SWENr1tIJ0C&printsec=frontcover&source=gbs_ge_summary_r&cad=0#v=onepage&q=sfu&f=false

"""

from __future__ import absolute_import, division, print_function
from astropy.table import Table
from sunpy.extern.six import  iteritems

from sunpy.sun import _constants as _con # pylint: disable=E0611

__all__ = ['get', 'find', 'print_all']

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
    <Constant name=u'Solar mass' value=1.9891e+30 uncertainty=5e+25 unit='kg' reference=u"Allen's Astrophysical Quantities 4th Ed.">
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
        data_rows.append([key, this_constant.name, this_constant.value, this_constant.uncertainty,
                          str(this_constant.unit), this_constant.reference])

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
