from __future__ import division

import numpy as np
from astropy import units as u
from astropy.coordinates import Longitude, Latitude, Angle
from sunpy.time import parse_time, julian_day

from sunpy.wcs import convert_hpc_hg, convert_hg_hpc
from sunpy.sun import constants, sun

__author__ = ["Jose Ivan Campos Rozo", "Stuart Mumford", "Jack Ireland"]
__all__ = ['diff_rot', 'rot_hpc']


@u.quantity_input(duration=u.s, latitude=u.degree)
def diff_rot(duration, latitude, rot_type='howard', frame_time='sidereal'):
    """
    This function computes the change in longitude over days in degrees.

    Parameters
    -----------
    duration : `~astropy.units.Quantity`
        Number of seconds to rotate over.
    latitude : `~astropy.units.Quantity`
        heliographic coordinate latitude in Degrees.
    rot_type : {'howard' | 'snodgrass' | 'allen'}
        howard : Use values for small magnetic features from Howard et al.
        snodgrass : Use Values from Snodgrass et. al
        allen : Use values from Allen, Astrophysical Quantities, and simpler equation.
    frame_time : {'sidereal' | 'synodic'}
        Choose 'type of day' time reference frame.

    Returns
    -------
    longitude_delta : `~astropy.units.Quantity`
        The change in longitude over days (units=degrees)

    Notes
    -----
    * IDL code equivalent: http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/diff_rot.pro
    * Howard rotation: http://adsabs.harvard.edu/abs/1990SoPh..130..295H
    * A review of rotation parameters (including Snodgrass values): http://link.springer.com/article/10.1023%2FA%3A1005226402796

    Examples
    --------
    Default rotation calculation over two days at 30 degrees latitude:

    >>> import numpy as np
    >>> import astropy.units as u
    >>> from sunpy.physics.differential_rotation import diff_rot
    >>> rotation = diff_rot(2 * u.day, 30 * u.deg)

    Default rotation over two days for a number of latitudes:

    >>> rotation = diff_rot(2 * u.day, np.linspace(-70, 70, 20) * u.deg)

    With rotation type 'allen':

    >>> rotation = diff_rot(2 * u.day, np.linspace(-70, 70, 20) * u.deg, 'allen')
    """

    latitude = latitude.to(u.deg)
    delta_seconds = duration.to(u.s).value
    delta_days = delta_seconds / 24.0 / 3600.0

    sin2l = (np.sin(latitude))**2
    sin4l = sin2l**2

    rot_params = {'howard': [2.894, -0.428, -0.370],
                  'snodgrass': [2.851, -0.343, -0.474]
                  }

    if rot_type not in ['howard', 'allen', 'snodgrass']:
        raise ValueError("""rot_type must equal one of
                        { 'howard' | 'allen' | 'snodgrass' }""")

    elif rot_type == 'allen':
        rotation_deg = delta_days * (14.44 - (3.0 * sin2l))

    else:
        A, B, C = rot_params[rot_type]

        # This is in micro-radians / sec
        rotation_rate = A + B * sin2l + C * sin4l
        rotation_deg = rotation_rate * 1e-6 * delta_seconds / np.deg2rad(1)

    if frame_time == 'synodic':
        rotation_deg -= 0.9856 * delta_days

    #return Longitude((np.round(rotation_deg, 4)), u.deg)
    return np.round(rotation_deg, 4) * u.deg


@u.quantity_input(x=u.arcsec, y=u.arcsec)
def rot_hpc(x, y, tstart, tend, frame_time='synodic', rot_type='howard', **kwargs):
    """Given a location on the Sun referred to using the Helioprojective
    Cartesian co-ordinate system (typically quoted in the units of arcseconds)
    use the solar rotation profile to find that location at some later or
    earlier time.  Note that this function assumes that the data was observed
    from the Earth or near Earth vicinity.  Specifically, data from SOHO and
    STEREO observatories are not supported.  Note also that the function does
    NOT use solar B0 and L0 values provided in source FITS files - these
    quantities are calculated.

    Parameters
    ----------
    x : `~astropy.units.Quantity`
        Helio-projective x-co-ordinate in arcseconds (can be an array).

    y : `~astropy.units.Quantity`
        Helio-projective y-co-ordinate in arcseconds (can be an array).

    tstart : `sunpy.time.time`
        date/time to which x and y are referred.

    tend : `sunpy.time.time`
    date/time at which x and y will be rotated to.

    rot_type : {'howard' | 'snodgrass' | 'allen'}
        | howard: Use values for small magnetic features from Howard et al.
        | snodgrass: Use Values from Snodgrass et. al
        | allen: Use values from Allen, Astrophysical Quantities, and simpler
          equation.

    frame_time : {'sidereal' | 'synodic'}
        Choose type of day time reference frame.

    Returns
    -------
    x : `~astropy.units.Quantity`
        Rotated helio-projective x-co-ordinate in arcseconds (can be an array).

    y : `~astropy.units.Quantity`
        Rotated helio-projective y-co-ordinate in arcseconds (can be an array).

    Examples
    --------
    >>> import astropy.units as u
    >>> from sunpy.physics.differential_rotation import rot_hpc
    >>> rot_hpc( -570 * u.arcsec, 120 * u.arcsec, '2010-09-10 12:34:56', '2010-09-10 13:34:56')
    (<Angle -562.9105822671319 arcsec>, <Angle 119.31920621992195 arcsec>)

    Notes
    -----
    SSWIDL code equivalent: http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/rot_xy.pro .
    The function rot_xy uses arcmin2hel.pro and hel2arcmin.pro to implement the
    same functionality as this function.  These two functions seem to perform
    inverse operations of each other to a high accuracy.  The corresponding
    equivalent functions here are convert_hpc_hg and convert_hg_hpc
    respectively. These two functions seem to perform inverse
    operations of each other to a high accuracy.  However, the values
    returned by arcmin2hel.pro are slightly different from those provided
    by convert_hpc_hg.  This leads to very slightly different results from
    rot_hpc compared to rot_xy.
    """

    # must have pairs of co-ordinates
    if np.array(x).shape != np.array(y).shape:
        raise ValueError('Input co-ordinates must have the same shape.')

    # Make sure we have enough time information to perform a solar differential
    # rotation
    # Start time
    dstart = parse_time(tstart)
    dend = parse_time(tend)
    interval = (dend - dstart).total_seconds() * u.s

    # Get the Sun's position from the vantage point at the start time
    vstart = kwargs.get("vstart", _calc_P_B0_SD(dstart))
    # Compute heliographic co-ordinates - returns (longitude, latitude). Points
    # off the limb are returned as nan
    longitude, latitude = convert_hpc_hg(x.to(u.arcsec).value,
                                         y.to(u.arcsec).value,
                                         b0_deg=vstart["b0"].to(u.deg).value,
                                         l0_deg=vstart["l0"].to(u.deg).value,
                                         dsun_meters=(constants.au * sun.sunearth_distance(t=dstart)).value,
                                         angle_units='arcsec')
    longitude = Longitude(longitude, u.deg)
    latitude = Angle(latitude, u.deg)
    # Compute the differential rotation
    drot = diff_rot(interval, latitude, frame_time=frame_time,
                    rot_type=rot_type)

    # Convert back to heliocentric cartesian in units of arcseconds
    vend = kwargs.get("vend", _calc_P_B0_SD(dend))

    # It appears that there is a difference in how the SSWIDL function
    # hel2arcmin and the sunpy function below performs this co-ordinate
    # transform.
    newx, newy = convert_hg_hpc(longitude.to(u.deg).value + drot.to(u.deg).value,
                                latitude.to(u.deg).value,
                                b0_deg=vend["b0"].to(u.deg).value,
                                l0_deg=vend["l0"].to(u.deg).value,
                                dsun_meters=(constants.au * sun.sunearth_distance(t=dend)).value,
                                occultation=False)
    newx = Angle(newx, u.arcsec)
    newy = Angle(newy, u.arcsec)
    return newx.to(u.arcsec), newy.to(u.arcsec)


def _calc_P_B0_SD(date):
    """
    To calculate the solar P, B0 angles and the semi-diameter as seen from
    Earth.  This function is assigned as being internal as these quantities
    should be calculated in a part of SunPy that can calculate these quantities
    accurately.

    Parameters
    -----------
    date : `sunpy.time.time`
        the time at which to calculate the solar P, B0 angles and the
        semi-diameter.

    Returns
    -------
    A dictionary with the following keys with the following meanings:

    p  -  Solar P (position angle of pole)  (degrees)
    b0 -  latitude of point at disk centre (degrees)
    sd -  semi-diameter of the solar disk in arcminutes

    Notes
    -----
    SSWIDL code equivalent:
        http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/pb0r.pro
    """
    # number of Julian days since 2415020.0
    de = julian_day(parse_time(date)) - 2415020.0

    # get the longitude of the sun etc.
    sun_position = _sun_pos(date)
    longmed = sun_position["longitude"].to(u.deg).value
    #ra = sun_position["ra"]
    #dec = sun_position["dec"]
    appl = sun_position["app_long"].to(u.deg).value
    oblt = sun_position["obliq"].to(u.deg).value

    # form the aberrated longitude
    Lambda = longmed - (20.50 / 3600.0)

    # form longitude of ascending node of sun's equator on ecliptic
    node = 73.6666660 + (50.250 / 3600.0) * ((de / 365.250) + 50.0)
    arg = Lambda - node

    # calculate P, the position angle of the pole
    p = np.rad2deg(
        np.arctan(-np.tan(np.deg2rad(oblt)) * np.cos(np.deg2rad(appl))) +
        np.arctan(-0.127220 * np.cos(np.deg2rad(arg))))

    # B0 the tilt of the axis...
    b = np.rad2deg(np.arcsin(0.12620 * np.sin(np.deg2rad(arg))))

    # ... and the semi-diameter
    # Form the mean anomalies of Venus(MV),Earth(ME),Mars(MM),Jupiter(MJ)
    # and the mean elongation of the Moon from the Sun(D).
    t = de / 36525.0
    mv = 212.60 + np.mod(58517.80 * t, 360.0)
    me = 358.4760 + np.mod(35999.04980 * t, 360.0)
    mm = 319.50 + np.mod(19139.860 * t, 360.0)
    mj = 225.30 + np.mod(3034.690 * t, 360.0)
    d = 350.70 + np.mod(445267.110 * t, 360.0)

    # Form the geocentric distance(r) and semi-diameter(sd)
    r = 1.0001410 - (0.0167480 - 0.00004180 * t) * np.cos(np.deg2rad(me)) \
        - 0.000140 * np.cos(np.deg2rad(2.0 * me)) \
        + 0.0000160 * np.cos(np.deg2rad(58.30 + 2.0 * mv - 2.0 * me)) \
        + 0.0000050 * np.cos(np.deg2rad(209.10 + mv - me)) \
        + 0.0000050 * np.cos(np.deg2rad(253.80 - 2.0 * mm + 2.0 * me)) \
        + 0.0000160 * np.cos(np.deg2rad(89.50 - mj + me)) \
        + 0.0000090 * np.cos(np.deg2rad(357.10 - 2.0 * mj + 2.0 * me)) \
        + 0.0000310 * np.cos(np.deg2rad(d))

    sd_const = constants.radius / constants.au
    sd = np.arcsin(sd_const / r) * 10800.0 / np.pi

    return {"p": Angle(p, u.deg),
            "b0": Angle(b, u.deg),
            "sd": Angle(sd.value, u.arcmin),
            "l0": Angle(0.0, u.deg)}


def _sun_pos(date):
    """
    Calculate solar ephemeris parameters.  Allows for planetary and lunar
    perturbations in the calculation of solar longitude at date and various
    other solar positional parameters. This routine is a truncated version of
    Newcomb's Sun and is designed to give apparent angular coordinates (T.E.D)
    to a precision of one second of time.  This function replicates the SSW/
    IDL function "sun_pos.pro".  This function is assigned to be
    internal at the moment as it should really be replaced by accurate
    ephemeris calculations in the part of SunPy that handles ephemeris.

    Parameters
    -----------
    date : `sunpy.time.time`
        Time at which the solar ephemeris parameters are calculated.  The
        input time can be in any acceptable time format.

    Returns
    -------
    A dictionary with the following keys with the following meanings:

    longitude  -  Longitude of sun for mean equinox of date (degs)
    ra         -  Apparent RA for true equinox of date (degs)
    dec        -  Apparent declination for true equinox of date (degs)
    app_long   -  Apparent longitude (degs)
    obliq      -  True obliquity (degs)

    Notes
    -----
    SSWIDL code equivalent:
        http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/sun_pos.pro

    Examples
    --------
    >>> from sunpy.physics.differential_rotation import _sun_pos
    >>> sp = _sun_pos('2013-03-27')
    """
    # Fractional Julian day with correct offset
    dd = julian_day(date) - 2415020.0

    # form time in Julian centuries from 1900.0
    t = dd / 36525.0

    # form sun's mean longitude
    l = (279.6966780 + np.mod(36000.7689250 * t, 360.00)) * 3600.0

    # allow for ellipticity of the orbit (equation of centre) using the Earth's
    # mean anomaly ME
    me = 358.4758440 + np.mod(35999.049750 * t, 360.0)
    ellcor = (6910.10 - 17.20 * t) * np.sin(np.deg2rad(me)) + \
    72.30 * np.sin(np.deg2rad(2.0 * me))
    l = l + ellcor

    # allow for the Venus perturbations using the mean anomaly of Venus MV
    mv = 212.603219 + np.mod(58517.8038750 * t, 360.0)
    vencorr = 4.80 * np.cos(np.deg2rad(299.10170 + mv - me)) + \
          5.50 * np.cos(np.deg2rad(148.31330 + 2.0 * mv - 2.0 * me)) + \
          2.50 * np.cos(np.deg2rad(315.94330 + 2.0 * mv - 3.0 * me)) + \
          1.60 * np.cos(np.deg2rad(345.25330 + 3.0 * mv - 4.0 * me)) + \
          1.00 * np.cos(np.deg2rad(318.150 + 3.0 * mv - 5.0 * me))
    l = l + vencorr

    # Allow for the Mars perturbations using the mean anomaly of Mars MM
    mm = 319.5294250 + np.mod(19139.858500 * t, 360.0)
    marscorr = 2.0 * np.cos(np.deg2rad(343.88830 - 2.0 * mm + 2.0 * me)) + \
            1.80 * np.cos(np.deg2rad(200.40170 - 2.0 * mm + me))
    l = l + marscorr

    # Allow for the Jupiter perturbations using the mean anomaly of Jupiter MJ
    mj = 225.3283280 + np.mod(3034.69202390 * t, 360.00)
    jupcorr = 7.20 * np.cos(np.deg2rad(179.53170 - mj + me)) + \
          2.60 * np.cos(np.deg2rad(263.21670 - mj)) + \
          2.70 * np.cos(np.deg2rad(87.14500 - 2.0 * mj + 2.0 * me)) + \
          1.60 * np.cos(np.deg2rad(109.49330 - 2.0 * mj + me))
    l = l + jupcorr

    # Allow for the Moons perturbations using the mean elongation of the Moon
    # from the Sun D
    d = 350.73768140 + np.mod(445267.114220 * t, 360.0)
    mooncorr = 6.50 * np.sin(np.deg2rad(d))
    l = l + mooncorr

    # Note the original code is
    # longterm  = + 6.4d0 * sin(( 231.19d0  +  20.20d0 * t )*!dtor)
    longterm = 6.40 * np.sin(np.deg2rad(231.190 + 20.20 * t))
    l = l + longterm
    l = np.mod(l + 2592000.0, 1296000.0)
    longmed = l / 3600.0

    # Allow for Aberration
    l = l - 20.5

    # Allow for Nutation using the longitude of the Moons mean node OMEGA
    omega = 259.1832750 - np.mod(1934.1420080 * t, 360.0)
    l = l - 17.20 * np.sin(np.deg2rad(omega))

    # Form the True Obliquity
    oblt = 23.4522940 - 0.01301250 * t + \
    (9.20 * np.cos(np.deg2rad(omega))) / 3600.0

    # Form Right Ascension and Declination
    l = l / 3600.0
    ra = np.rad2deg(np.arctan2(np.sin(np.deg2rad(l)) * \
                        np.cos(np.deg2rad(oblt)), np.cos(np.deg2rad(l))))

    if isinstance(ra, np.ndarray):
        ra[ra < 0.0] += 360.0
    elif ra < 0.0:
        ra = ra + 360.0

    dec = np.rad2deg(np.arcsin(np.sin(np.deg2rad(l)) *
                               np.sin(np.deg2rad(oblt))))

    # convert the internal variables to those listed in the top of the
    # comment section in this code and in the original IDL code.  Quantities
    # are assigned following the advice in Astropy "Working with Angles"
    return {"longitude": Longitude(longmed, u.deg),
            "ra": Longitude(ra, u.deg),
            "dec": Latitude(dec, u.deg),
            "app_long": Longitude(l, u.deg),
            "obliq": Angle(oblt, u.deg)}

