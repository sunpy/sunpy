"""
Fundamental Solar Physical Constants
------------------------------------
These constants are taken from various sources. The structure of this module is heavily 
based on if not directly copied from the SciPy constants module but contains Solar 
Physical constants. All units are in SI (mks) unless otherwise specified.

Object
------
physical_constants : dict
    A dictionary containing physical constants. Keys are the names
    of physical constants, values are tuples (value, units, precision). The dictionary
    contains the following solar physical constants:
        
    GM: 
        The gravitational constant multiplied by the mass of the Sun.
    absolute magnitude: 
        The absolute visual magnitude of the Sun. It is the measure of the Sun's intrinsic 
        brightness. It is the apparent magnitude the Sun would have if it were 32.6 light
        years (10 parsecs) away from Earth.
    average density:
         The average density of the Sun in SI.
    average_angular_size: 
        The average angular size of the Sun as seen from Earth in arcseconds.
    center density:
        The density at the center of the Sun in SI.
    center temperature:
        The temperature at the center of the Sun in Kelvin.
    diameter: 
        The diameters of the Sun at the equator in meters.
    effective temperature:
        The effective black-body temperature of the Sun in Kelvin.  
    ellipticity: 
        The ellipticity of the Sun.
    escape velocity: 
        The velocity which an object needs to escape from the gravitational pull of the Sun.
    luminosity:
        The luminosity of the Sun in Joules per second.
    mass:
        The mass of the Sun in kg.
    mass conversion rate:
        The rate at which the Sun converts mass to energy.
    mean energy production:
        The mean rate at which the Sun produces energy in Joules per second
    mean intensity:
        The mean intensity of the Sun.
    metallicity: 
        The metallicity of the Sun
    radius:
        The radius of the Sun at the equator in meters.
    solar flux unit:
        The definition of a solar flux unit.
    sunspot cycle:
        The average duration of the solar activity cycle.
    surface area:
        The surface area of the Sun in meters squared. 
    surface gravity:
        The gravitational acceleration at the surface of the Sun as measured at the equator.
    visual magnitude:
       A measure of the Sun's brightness as seen by an observer on Earth without the
       presence of the atmosphere.
    volume:
        The volume of the Sun in meters cubed.

Attributes
----------
A number of variables from physical_constants are made available for convenience as 
attributes. They are equatorial_radius, radius (same as equatorial radius), equatorial_diameter, volume, surface_area, average_density, center_density, 
equatorial_surface_gravity, mean_intensity, effective_temperature, center_temperature, 
luminosity, absolute_magnitude, visual_magnitude, mass_conversion_rate, 
mean_energy_production, escape_velocity, ellipticity, GM, average_angular_size, sfu.

Source
------
Constants are imported from Review of Particle Physics 2010 (page 102), 
and NASA's Sun Fact Sheet as well as other sources.

Websites
--------
http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html

.. todo:: Need better sources as well as error values.

.. todo:: Create a cheat sheet function which prints out key solar values.

"""
#TODO: Need better sources as well as error values.
#TODO: Create a cheat sheet function which prints out key solar values.

from __future__ import absolute_import

import scipy.constants as _cd
from sunpy.sun import _si as _con # pylint: disable=E0611

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


def print_all(key = None):
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
