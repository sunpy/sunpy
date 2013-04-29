"""
Position of the Sun.
"""
from __future__ import absolute_import

import numpy as np
#pylint: disable=E1101,E1121

__all__ = ['position']

__authors__ = ["Jack Ireland"]
__email__ = "jack.ireland@nasa.gov"

"""
Questions (now an issue: https://github.com/sunpy/sunpy/issues/394)
---------
1. should the output always be an ndarray?
2. Should we be able to return position in other coordinate systems as well?
   If so, then this function should be to something like radec, and pos() should
   accept "radec" as input and call the current function.
3. Where do all of the "magic" numbers come from? If possible document them, and
   they may be useful elsewhere, consider putting them in sunpy.constants (could
   any of them already exist there or in scipy constants?)
"""

def position(date, radian=False):
    """
    Routine to calculate the right ascension (RA) and declination (dec) of 
    the Sun.
    
    Based in Solarsoft/IDL routine SUNPOS (Sept 1997)
    
    Parameters
    ----------
    date : scalar, numpy.ndarray ???
        Julian dates as input
    radian : boolean
        if True, return the results in radians

    Attributes
    ----------
    none

    Examples
    --------
    >>> from sunpy.sun import pos
    >>> pos.pos(0.0)
    >>> (array([ 241.65373812]), array([-21.68562104]), array([ 243.82420808]), array([ 24.31513199]))'
    >>> pos.pos(2455813.024259259, radian=True)
    >>> (array([ 2.90902792]), array([ 0.09958257]), array([ 2.88897706]), array([ 0.40905761]))
    
    See Also
    --------
    numpy.ndarray
    
    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
    """
    
    # Make sure the input is a ndarray
    if np.isscalar(date):
        date = date + np.zeros([1])

    # Degrees to radians
    deg2rad = np.pi / 180.0
    
    # Time in Julian centuries from 1900.0
    time = (date - 2415020.0) / 36525.0

    # Sun's mean longitude
    longitude = (279.696678 + np.mod(36000.768925 * time, 360.0)) * 3600.0
    
    # Allow for ellipticity of the orbit (equation of centre)
    # using the Earth's mean anomaly ME
    me = 358.4758440 + np.mod(35999.0497500 * time, 360.0)
    ellcor  = ((6910.10 - 17.20 * time) * np.sin(me * deg2rad) + 
               72.30 * np.sin(2.00 * me * deg2rad))
    longitude = longitude + ellcor

    # Allow for the Venus perturbations using the mean anomaly of Venus MV
    mv = 212.6032190 + np.mod(58517.8038750 * time, 360.0)
    vencorr = (4.80 * np.cos((299.10170 + mv - me) * deg2rad) + 
               5.50 * np.cos((148.31330 + 2.00 * mv - 2.00 * me ) * deg2rad) + 
               2.50 * np.cos((315.94330 + 2.00 * mv - 3.00 * me ) * deg2rad) + 
               1.60 * np.cos((345.25330 + 3.00 * mv - 4.00 * me ) * deg2rad) + 
               1.00 * np.cos((318.150   + 3.00 * mv - 5.00 * me ) * deg2rad))
    longitude = longitude + vencorr

    # Allow for the Mars perturbations using the mean anomaly of Mars MM
    mm = 319.5294250 + np.mod( 19139.8585000 * time,  360.00)
    marscorr = (2.00 * np.cos((343.88830 - 2.00 * mm + 2.00 * me) * deg2rad ) + 
                1.80 * np.cos((200.40170 - 2.00 * mm + me) * deg2rad))
    longitude = longitude + marscorr

    # Allow for the Jupiter perturbations using the mean anomaly of
    # Jupiter MJ
    mj = 225.3283280 + np.mod( 3034.69202390 * time, 360.00 )
    jupcorr = (7.20 * np.cos(( 179.53170 - mj + me ) * deg2rad) +
               2.60 * np.cos((263.21670 - mj) * deg2rad) +
               2.70 * np.cos(( 87.14500 - 2.00 * mj + 2.00 * me) * deg2rad) +
               1.60 * np.cos((109.49330 - 2.00 * mj + me) * deg2rad))
    longitude = longitude + jupcorr
    
    # Allow for the Moons perturbations using the mean elongation of
    # the Moon from the Sun
    moon_me = 350.73768140 + np.mod(445267.114220 * time, 360.00)
    mooncorr = 6.50 * np.sin(moon_me * deg2rad)
    longitude = longitude + mooncorr

    # Allow for long period terms
    longterm = + 6.40 * np.sin((231.190 + 20.200 * time) * deg2rad)
    longitude = longitude + longterm
    longitude = np.mod(longitude + 2592000.00, 1296000.00)
    longmed = longitude / 3600.00

    # Allow for Aberration ??????
    longitude = longitude - 20.50

    # Allow for Nutation using the longitude of the Moons mean node OMEGA
    omega = 259.1832750 - np.mod(1934.1420080 * time, 360.00)
    longitude = longitude - 17.20 * np.sin(omega * deg2rad)

    # Form the True Obliquity
    obliquity = (23.4522940 - 0.01301250 * 
                 time + (9.20 * np.cos(omega * deg2rad)) / 3600.00)

    # Form Right Ascension and Declination
    longitude = longitude / 3600.00
    ra = np.arctan2(np.sin(longitude * deg2rad) * np.cos(obliquity * deg2rad), 
                    np.cos(longitude * deg2rad))

    # ?
    neg = np.where(ra < 0.0)
    if len(neg) > 0:
        ra[neg] = ra[neg] + 2.0 * np.pi

    dec = np.arcsin(np.sin(longitude * deg2rad) * np.sin(obliquity * deg2rad))

    if radian:
        obliquity = obliquity * deg2rad 
        longmed = longmed * deg2rad
    else:
        ra = ra / deg2rad
        dec = dec / deg2rad

    return ra, dec, longmed, obliquity