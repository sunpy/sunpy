"""
Fundamental Solar Physical Constants
------------------------------------
These constants are taken from various sources. The structure of this module is heavily
based on if not directly copied from the SciPy constants module but contains Solar
Physical constants.

Object
------
    physical_constants : dict
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
A number of variables from physical_constants are made available for convenience as
attributes.

Websites
--------
| http://books.google.com/books?id=4SWENr1tIJ0C&printsec=frontcover&source=gbs_ge_summary_r&cad=0#v=onepage&q=sfu&f=false

"""

from __future__ import absolute_import

from sunpy.sun import _constants as _con # pylint: disable=E0611

physical_constants = _con.physical_constants


def constant(key):
    """
    The constant in physical_constants index by key

    Parameters
    ----------
    key : Python string or unicode
        Key in dictionary in `physical_constants`

    Returns
    -------
    constant : constant
        constant in `physical_constants` corresponding to `key`

    See Also
    --------
    _constants : Contains the description of `physical_constants`, which, as a
        dictionary literal object, does not itself possess a docstring.

    Examples
    --------
    >>> from sunpy.sun import constants
    >>> constants.constant('mass')
    <Constant name=u'Solar mass' value=1.9891e+30 error=5e+25 units='kg' reference=u"Allen's Astrophysical Quantities 4th Ed.">
    """
    return physical_constants[key]


def value(key):
    """
    Value in physical_constants indexed by key

    Parameters
    ----------
    key : Python string or unicode
        Key in dictionary `physical_constants`

    Returns
    -------
    value : float
        Value in `physical_constants` corresponding to `key`

    See Also
    --------
    _constants : Contains the description of `physical_constants`, which, as a
        dictionary literal object, does not itself possess a docstring.

    Examples
    --------
    >>> from sunpy.sun import constants
    >>> constants.uncertainty('mass')
    5e+25

    """
    return constant(key).value


def unit(key):
    """
    Unit in physical_constants indexed by key

    Parameters
    ----------
    key : Python string or unicode
        Key in dictionary `physical_constants`

    Returns
    -------
    unit : Python string
        Unit in `physical_constants` corresponding to `key`

    See Also
    --------
    _constants : Contains the description of `physical_constants`, which, as a
        dictionary literal object, does not itself possess a docstring.

    Examples
    --------
    >>> from sunpy.sun import constants
    >>> constants.uncertainty('mass')
    5e+25

    """
    return constant(key).unit


def uncertainty(key):
    """
    Relative uncertainty in physical_constants indexed by key

    Parameters
    ----------
    key : Python string or unicode
        Key in dictionary `physical_constants`

    Returns
    -------
    prec : float
        Relative uncertainty in `physical_constants` corresponding to `key`

    See Also
    --------
    _constants : Contains the description of `physical_constants`, which, as a
        dictionary literal object, does not itself possess a docstring.

    Examples
    --------
    >>> from sunpy.sun import constants
    >>> constants.uncertainty('mass')
    5e+25

    """
    return constant(key).uncertainty


def find(sub=None, disp=False):
    """
    Return list of physical_constants keys containing a given string

    Parameters
    ----------
    sub : str, unicode
        Sub-string to search keys for.  By default, return all keys.
    disp : bool
        If True, print the keys that are found, and return None.
        Otherwise, return the list of keys without printing anything.

    Returns
    -------
    keys : None or list
        If `disp` is False, the list of keys is returned. Otherwise, None
        is returned.

    See Also
    --------
    _constants : Contains the description of `physical_constants`, which, as a
        dictionary literal object, does not itself possess a docstring.

    """
    if sub is None:
        result = physical_constants.keys()
    else:
        result = [key for key in physical_constants \
                 if sub.lower() in key.lower()]

    result.sort()
    if disp:
        for key in result:
            print key
        return
    else:
        return result


def print_all(key=None):
    """
    Prints out the complete list of physical_constants to the screen or
    one single value

    Parameters
    ----------
    key : Python string or unicode
        Key in dictionary `physical_constants`

    Returns
    -------
    None

    See Also
    --------
    _constants : Contains the description of `physical_constants`, which, as a
        dictionary literal object, does not itself possess a docstring.

    """
    column_width = [25, 20, 20, 20]
    table_width = (column_width[0] + column_width[1] + column_width[2]
                   + column_width[3])
    format_string = ('{0:<' + str(column_width[0]) + '}' + '{1:>' +
                    str(column_width[1]) + '}' + '{2:>' + str(column_width[2])
                    + '}' + '{3:>' + str(column_width[3]) + '}')
    print(format_string.format('Name', 'Value', 'Units', 'Error'))
    print(('{:-^' + str(table_width) + '}').format(''))

    if key is None:
        for key in physical_constants:
            print(format_string.format(key, str(value(key)), unit(key),
                                       str(uncertainty(key))))
    else:
        print(format_string.format(key, str(value(key)), unit(key),
                                   str(uncertainty(key))))

# Spectral class is not included in physical constants since it is not a number
spectral_classification = 'G2V'

au = astronomical_unit = constant('mean distance')

# The following variables from _constants are brought out by making them
# accessible through a call such as sun.volume
mass = constant('mass')
equatorial_radius = radius = constant('radius')
volume = constant('volume')
surface_area = constant('surface area')
average_density = density = constant('average density')
equatorial_surface_gravity = surface_gravity = constant('surface gravity')
effective_temperature = constant('effective temperature')
luminosity = constant('luminosity')
mass_conversion_rate = constant('mass conversion rate')
escape_velocity = constant('escape velocity')

sfu = constant('solar flux unit')

# Observable parameters
average_angular_size = constant('average angular size')
