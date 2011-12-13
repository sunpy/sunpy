"""
Fundamental Solar Physical Constants
------------------------------
These constants are taken from various sources. This module is heavily based on
(if not directly copied from) SciPy constants module.

Object
------
physical_constants : dict
    A dictionary containing physical constants. Keys are the names
    of physical constants, values are tuples (value, units, precision).

Functions
---------
value(key):
    Returns the value of the physical constant(key).
unit(key):
    Returns the units of the physical constant(key).
precision(key):
    Returns the relative precision of the physical constant(key).
find(sub):
    Prints or returns list of keys containing the string sub,
    default is all.

Source
------
Constants are imported from Review of Particle Physics 2010 (page 102), 
and NASA's Sun Fact Sheet as well as other sources.

Websites
-------
http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html

TODO: 
----------
Need better sources as well as error values.
Create a cheat sheet function which prints out key solar values
"""
from __future__ import absolute_import

import scipy.constants as _cd
from . import _si as _con

physical_constants = _con.physical_constants

au = astronomical_unit = _cd.au

# The following functions (value, precision, unit, find) are copied directly 
# from SciPy constants.
def value(key) :
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
    >>> constants.precision('mass')
        1.9884e30

    """
    return physical_constants[key][0]

def unit(key) :
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
    >>> constants.precision('mass')
    'kg'

    """
    return physical_constants[key][1]

def precision(key) :
    """
    Relative precision in physical_constants indexed by key

    Parameters
    ----------
    key : Python string or unicode
        Key in dictionary `physical_constants`

    Returns
    -------
    prec : float
        Relative precision in `physical_constants` corresponding to `key`

    See Also
    --------
    _constants : Contains the description of `physical_constants`, which, as a
        dictionary literal object, does not itself possess a docstring.

    Examples
    --------
    >>> from sunpy.sun import constants
    >>> constants.precision('mass')
    

    """
    return physical_constants[key][2] / physical_constants[key][0]

def find(sub=None, disp=False):
    """
    Return list of physical_constant keys containing a given string

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


def print_all(key = None):
    """
    Prints out the complete list of physical_constant to the screen or
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
    print(format_string.format('Name', 'Value', 'Units', 'Precision'))
    print(('{:-^' + str(table_width) + '}').format(''))

    if key is None:
        for key in physical_constants:
            print(format_string.format(key, str(value(key)), unit(key), 
                                       str(precision(key))))
    else: 
            print(format_string.format(key, str(value(key)), unit(key), 
                                       str(precision(key))))

# Spectral class is not included in physical constants since it is not a number
spectral_classification = 'G2V'

# The following variables from _constants are brought out by making them 
# accessible through a call such as sun.volume
equatorial_radius = radius = value('radius')
equatorial_diameter = value('diameter')
volume = value('volume')
surface_area = value('surface area')
average_density = density = value('average density')
center_density = value('center density')
equatorial_surface_gravity = surface_gravity = value('surface gravity')
mean_intensity = intensity = value('mean intensity')
effective_temperature = value('effective temperature')
center_temperature = value('center temperature')
luminosity = value('luminosity')
absolute_magnitude = value('absolute magnitude')
visual_magnitude = value('visual magnitude')
mass_conversion_rate = value('mass conversion rate')
mean_energy_production = value('mean energy production')
escape_velocity = value('escape velocity')
ellipticity = value('ellipticity')
GM = value('GM')

sfu = value('solar flux unit')

# Observable parameters
average_angular_size = value('average_angular_size')
