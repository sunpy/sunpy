from __future__ import absolute_import
import numpy as np

"""
Position of the Sun in the sky as seen from Earth.
"""

#pylint: disable=E1101,E1121
__authors__ = ["Jack Ireland"]
__email__ = "jack.ireland@nasa.gov"



"""
Questions
---------
1. should the output always be an ndarray?
"""


def pos(jd, radian = False):
    """
    pos(jd, radian = False)
    
    Routine to calculate the right ascension (RA) and declination (dec) of the Sun.
    
    Parameters
    ----------
    jd : scalar, numpy.ndarray
        Julian dates as input
    radian : boolean
        if True, return the results in radians

    Attributes
    ----------
    none

    Examples
    --------
    >>> sunpy.sun.pos(0.0)
    >>> (array([ 241.65373812]), array([-21.68562104]), array([ 243.82420808]), array([ 24.31513199]))'
    >>> sunpy.sun.pos(2455813.024259259,radian=True)
    >>> (array([ 2.90902792]), array([ 0.09958257]), array([ 2.88897706]), array([ 0.40905761]))
    
    See Also
    --------
    numpy.ndarray parent class
    
    References
    ----------
    | http://docs.scipy.org/doc/numpy/reference/arrays.classes.html
    | Original Solarsoft/IDL documentation follows:
    ;+
    ; NAME:
    ;       SUNPOS
    ; PURPOSE:
    ;       To compute the RA and Dec of the Sun at a given date.
    ;
    ; CALLING SEQUENCE: IDL
    ;       SUNPOS, jd, ra, dec, [elong, obliquity, /RADIAN ]
    ; INPUTS: IDL
    ;       jd    - The Julian date of the day (and time), scalar or vector
    ;               usually double precision
    ; OUTPUTS: IDL
    ;       ra    - The right ascension of the sun at that date in DEGREES
    ;               double precision, same number of elements as jd
    ;       dec   - The declination of the sun at that date in DEGREES
    ; OPTIONAL OUTPUTS: IDL
    ;       elong - Ecliptic longitude of the sun at that date in DEGREES.
    ;       obliquity - the obliquity of the ecliptic, in DEGREES
    ;
    ; OPTIONAL INPUT KEYWORD: IDL
    ;       /RADIAN - If this keyword is set and non-zero, then all output variables 
    ;               are given in Radians rather than Degrees
    ;
    ; NOTES:
    ;       Patrick Wallace (Rutherford Appleton Laboratory, UK) has tested the
    ;       accuracy of a C adaptation of the sunpos.pro code and found the 
    ;       following results.   From 1900-2100 SUNPOS  gave 7.3 arcsec maximum 
    ;       error, 2.6 arcsec RMS.  Over the shorter interval 1950-2050 the figures
    ;       were 6.4 arcsec max, 2.2 arcsec RMS.  
    ;
    ;       The returned RA and Dec are in the given date's equinox.
    ;
    ;       Procedure was extensively revised in May 1996, and the new calling
    ;       sequence is incompatible with the old one.
    ; METHOD:
    ;       Uses a truncated version of Newcomb's Sun.    Adapted from the IDL
    ;       routine SUN_POS by CD Pike, which was adapted from a FORTRAN routine
    ;       by B. Emerson (RGO).
    ; EXAMPLE:
    ;       (1) Find the apparent RA and Dec of the Sun on May 1, 1982
    ;       
    ;       IDL> jdcnv, 1982, 5, 1,0 ,jd      ;Find Julian date jd = 2445090.5   
    ;       IDL> sunpos, jd, ra, dec
    ;       IDL> print,adstring(ra,dec,2)
    ;                02 31 32.61  +14 54 34.9
    ;
    ;       The Astronomical Almanac gives 02 31 32.58 +14 54 34.9 so the error
    ;               in SUNPOS for this case is < 0.5 (double quote).      
    ;
    ;       (2) Find the apparent RA and Dec of the Sun for every day in 1997
    ;
    ;       IDL> jdcnv, 1997,1,1,0, jd                ;Julian date on Jan 1, 1997
    ;       IDL> sunpos, jd+ dindgen(365), ra, dec    ;RA and Dec for each day 
    ;
    ; MODIFICATION HISTORY:
    ;       Written by Michael R. Greason, STX, 28 October 1988.
    ;       Accept vector arguments, W. Landsman     April,1989
    ;       Eliminated negative right ascensions.  MRG, Hughes STX, 6 May 1992.
    ;       Rewritten using the 1993 Almanac.  Keywords added.  MRG, HSTX, 
    ;               10 February 1994.
    ;       Major rewrite, improved accuracy, always return values in degrees
    ;       W. Landsman  May, 1996 
    ;       Added /RADIAN keyword,    W. Landsman       August, 1997
    ;       Converted to IDL V5.0   W. Landsman   September 1997
    """
    
    # make sure the input is a ndarray
    if np.isscalar(jd):
        jd = jd + np.zeros([1])
    # degrees to radians
    dtor = np.pi/180.0
    # form time in Julian centuries from 1900.0
    t = (jd - 2415020.0)/36525.0

    #  form sun's mean longitude
    l = (279.696678 + np.mod( 36000.768925*t, 360.0 ))*3600.0
    
    #  allow for ellipticity of the orbit (equation of centre)
    #  using the Earth's mean anomaly ME
    me = 358.4758440 + np.mod(35999.0497500*t, 360.0)
    ellcor  = (6910.10 - 17.20*t)*np.sin(me*dtor) + 72.30*np.sin(2.00*me*dtor)
    l = l + ellcor

    # allow for the Venus perturbations using the mean anomaly of Venus MV
    mv = 212.6032190 + np.mod(58517.8038750*t, 360.0)
    vencorr = 4.80 * np.cos((299.10170 + mv - me)*dtor) + \
              5.50 * np.cos((148.31330 +  2.00 * mv  -  2.00 * me )*dtor) + \
              2.50 * np.cos((315.94330 +  2.00 * mv  -  3.00 * me )*dtor) + \
              1.60 * np.cos((345.25330 +  3.00 * mv  -  4.00 * me )*dtor) + \
              1.00 * np.cos((318.150   +  3.00 * mv  -  5.00 * me )*dtor)
    l = l + vencorr

    #  Allow for the Mars perturbations using the mean anomaly of Mars MM
    mm = 319.5294250  +  np.mod( 19139.8585000 * t,  360.00 )
    marscorr = 2.00 * np.cos((343.88830 -  2.00 * mm  +  2.00 * me)*dtor ) + \
               1.80 * np.cos((200.40170 -  2.00 * mm  + me) * dtor)
    l = l + marscorr

    # Allow for the Jupiter perturbations using the mean anomaly of
    # Jupiter MJ
    mj = 225.3283280  +  np.mod( 3034.69202390 * t,  360.00 )
    jupcorr = 7.20 * np.cos(( 179.53170 - mj + me )*dtor) + \
              2.60 * np.cos((263.21670  -  mj ) *dtor) + \
              2.70 * np.cos(( 87.14500  -  2.00 * mj  +  2.00 * me ) *dtor) + \
              1.60 * np.cos((109.49330  -  2.00 * mj  +  me ) *dtor)
    l = l + jupcorr
    
    # Allow for the Moons perturbations using the mean elongation of
    # the Moon from the Sun D
    d = 350.73768140  + np.mod( 445267.114220 * t,  360.00 )
    mooncorr  = 6.50 * np.sin(d*dtor)
    l = l + mooncorr

    # Allow for long period terms
    longterm  = + 6.40 * np.sin(( 231.190  +  20.200 * t )*dtor)
    l  =    l + longterm
    l  =  np.mod( l + 2592000.00, 1296000.00 )
    longmed = l/3600.00

    # Allow for Aberration
    l  =  l - 20.50

    # Allow for Nutation using the longitude of the Moons mean node OMEGA
    omega = 259.1832750 - np.mod( 1934.1420080 * t, 360.00 )
    l  =  l - 17.20 * np.sin(omega*dtor)

    # Form the True Obliquity
    oblt  = 23.4522940 - 0.01301250*t + (9.20*np.cos(omega*dtor))/3600.00

    # Form Right Ascension and Declination
    l = l/3600.00
    ra  = np.arctan2( np.sin(l*dtor) * np.cos(oblt*dtor) , np.cos(l*dtor) )

    #
    neg = np.where(ra<0.0)
    Nneg = len(neg)
    if Nneg > 0:
        ra[neg] = ra[neg] + 2.0*np.pi

    dec = np.arcsin( np.sin(l*dtor) * np.sin(oblt*dtor) )

    if radian:
        oblt = np.deg2rad(oblt) 
        longmed = np.deg2rad(longmed)
    else:
        ra = ra/dtor
        dec = dec/dtor

    return ra, dec, longmed, oblt