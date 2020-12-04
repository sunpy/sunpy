"""
Ephemeris calculations using SunPy coordinate frames
"""
import numpy as np

import astropy.units as u
from astropy.constants import c as speed_of_light
from astropy.coordinates import (
    ICRS,
    HeliocentricEclipticIAU76,
    SkyCoord,
    get_body_barycentric,
    get_body_barycentric_posvel,
)
from astropy.coordinates.representation import (
    CartesianDifferential,
    CartesianRepresentation,
    SphericalRepresentation,
)

from sunpy import log
from sunpy.time import parse_time
from sunpy.time.time import _variables_for_parse_time_docstring
from sunpy.util.decorators import add_common_docstring
from .frames import HeliographicStonyhurst

__author__ = "Albert Y. Shih"
__email__ = "ayshih@gmail.com"

__all__ = ['get_body_heliographic_stonyhurst', 'get_earth',
           'get_horizons_coord']


@add_common_docstring(**_variables_for_parse_time_docstring())
def get_body_heliographic_stonyhurst(body, time='now', observer=None, *, include_velocity=False):
    """
    Return a `~sunpy.coordinates.frames.HeliographicStonyhurst` frame for the location of a
    solar-system body at a specified time.  The location can be corrected for light travel time
    to an observer.

    Parameters
    ----------
    body : `str`
        The solar-system body for which to calculate positions
    time : {parse_time_types}
        Time to use in a parse_time-compatible format
    observer : `~astropy.coordinates.SkyCoord`
        If None, the returned coordinate is the instantaneous or "true" location.
        If not None, the returned coordinate is the astrometric location (i.e., accounts for light
        travel time to the specified observer)

    Keyword Arguments
    -----------------
    include_velocity : `bool`
        If True, include the body's velocity in the output coordinate.  Defaults to False.

    Returns
    -------
    out : `~sunpy.coordinates.frames.HeliographicStonyhurst`
        Location of the solar-system body in the `~sunpy.coordinates.HeliographicStonyhurst` frame

    Notes
    -----
    There is no correction for aberration due to observer motion.  For a body close to the Sun in
    angular direction relative to the observer, the correction can be negligible because the
    apparent location of the body will shift in tandem with the Sun.

    Examples
    --------
    >>> from sunpy.coordinates.ephemeris import get_body_heliographic_stonyhurst

    Obtain the location of Venus

    >>> get_body_heliographic_stonyhurst('venus', '2012-06-06 04:07:29')
    <HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
        (0.07349535, 0.05223575, 0.72605496)>

    Obtain the location of Venus as seen from Earth when adjusted for light travel time

    >>> earth = get_body_heliographic_stonyhurst('earth', '2012-06-06 04:07:29')
    >>> get_body_heliographic_stonyhurst('venus', '2012-06-06 04:07:29', observer=earth)
    INFO: Apparent body location accounts for 144.07 seconds of light travel time [sunpy.coordinates.ephemeris]
    <HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
        (0.07084926, 0.0520573, 0.72605477)>

    Obtain the location and velocity of Mars

    >>> mars = get_body_heliographic_stonyhurst('mars', '2001-02-03', include_velocity=True)
    >>> mars
    <HeliographicStonyhurst Coordinate (obstime=2001-02-03T00:00:00.000): (lon, lat, radius) in (deg, deg, AU)
        (63.03105777, -5.20656151, 1.6251161)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (-0.02323686, 0.00073376, -1.4798387)>

    Transform that same location and velocity of Mars to a different frame using
    `~astropy.coordinates.SkyCoord`.

    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import Helioprojective
    >>> SkyCoord(mars).transform_to(Helioprojective(observer=earth))
    <SkyCoord (Helioprojective: obstime=2001-02-03T00:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000): (lon, lat, radius) in (deg, deg, AU)
        (7.835757e-15, -0.00766698, 1.01475668)>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
        (-298029.94625805, -21753.50941181, 1.40010091)
     (d_Tx, d_Ty, d_distance) in (arcsec / s, arcsec / s, km / s)
        (-0.01652981, -0.00059216, -15.14320414)>
    """
    obstime = parse_time(time)

    if observer is None:
        # If there is no observer, there is not adjustment for light travel time
        emitted_time = obstime
    else:
        observer_icrs = SkyCoord(observer).icrs.cartesian

        # This implementation is modeled after Astropy's `_get_apparent_body_position`
        light_travel_time = 0.*u.s
        emitted_time = obstime
        delta_light_travel_time = 1.*u.s  # placeholder value
        while np.any(np.fabs(delta_light_travel_time) > 1.0e-8*u.s):
            body_icrs = get_body_barycentric(body, emitted_time)
            distance = (body_icrs - observer_icrs).norm()
            delta_light_travel_time = light_travel_time - distance / speed_of_light
            light_travel_time = distance / speed_of_light
            emitted_time = obstime - light_travel_time

        if light_travel_time.isscalar:
            ltt_string = f"{light_travel_time.to_value('s'):.2f}"
        else:
            ltt_string = f"{light_travel_time.to_value('s')}"
        log.info(f"Apparent body location accounts for {ltt_string} seconds of light travel time")

    if include_velocity:
        pos, vel = get_body_barycentric_posvel(body, emitted_time)
        body_icrs = pos.with_differentials(vel.represent_as(CartesianDifferential))
    else:
        body_icrs = get_body_barycentric(body, emitted_time)

    body_hgs = ICRS(body_icrs).transform_to(HeliographicStonyhurst(obstime=obstime))

    return body_hgs


@add_common_docstring(**_variables_for_parse_time_docstring())
def get_earth(time='now', *, include_velocity=False):
    """
    Return a `~astropy.coordinates.SkyCoord` for the location of the Earth at a specified time in
    the `~sunpy.coordinates.frames.HeliographicStonyhurst` frame.  The longitude will be zero by
    definition.

    Parameters
    ----------
    time : {parse_time_types}
        Time to use in a parse_time-compatible format

    Keyword Arguments
    -----------------
    include_velocity : `bool`
        If True, include the Earth's velocity in the output coordinate. Defaults to False.

    Returns
    -------
    out : `~astropy.coordinates.SkyCoord`
        Location of the Earth in the `~sunpy.coordinates.frames.HeliographicStonyhurst` frame

    Notes
    -----
    The Earth's velocity in the output coordinate will invariably be negligible in the longitude
    direction because the `~sunpy.coordinates.frames.HeliographicStonyhurst` frame rotates in time
    such that the plane of zero longitude (the XZ-plane) tracks Earth.

    Examples
    --------
    >>> from sunpy.coordinates.ephemeris import get_earth
    >>> get_earth('2001-02-03 04:05:06')
    <SkyCoord (HeliographicStonyhurst: obstime=2001-02-03T04:05:06.000): (lon, lat, radius) in (deg, deg, AU)
        (0., -6.18656962, 0.98567647)>
    >>> get_earth('2001-02-03 04:05:06', include_velocity=True)
    <SkyCoord (HeliographicStonyhurst: obstime=2001-02-03T04:05:06.000): (lon, lat, radius) in (deg, deg, AU)
        (0., -6.18656962, 0.98567647)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (6.42643739e-11, -0.00279484, 0.24968506)>
    >>> get_earth('2001-02-03 04:05:06', include_velocity=True).transform_to('heliocentricinertial')
    <SkyCoord (HeliocentricInertial: obstime=2001-02-03T04:05:06.000): (lon, lat, distance) in (deg, deg, AU)
        (58.41594489, -6.18656962, 0.98567647)
     (d_lon, d_lat, d_distance) in (arcsec / s, arcsec / s, km / s)
        (0.0424104, -0.00279484, 0.2496851)>
    """
    earth = get_body_heliographic_stonyhurst('earth', time=time, include_velocity=include_velocity)

    # Explicitly set the longitude to 0
    earth_repr = SphericalRepresentation(0*u.deg, earth.lat, earth.radius)

    # Modify the representation in the frame while preserving all differentials (e.g., velocity)
    earth = earth.realize_frame(earth_repr.with_differentials(earth.spherical.differentials))

    return SkyCoord(earth)


@add_common_docstring(**_variables_for_parse_time_docstring())
def get_horizons_coord(body, time='now', id_type='majorbody', *, include_velocity=False):
    """
    Queries JPL HORIZONS and returns a `~astropy.coordinates.SkyCoord` for the location of a
    solar-system body at a specified time.  This location is the instantaneous or "true" location,
    and is not corrected for light travel time or observer motion.

    .. note::
        This function requires the Astroquery package to be installed and
        requires an Internet connection.

    Parameters
    ----------
    body : `str`
        The solar-system body for which to calculate positions.  One can also use the search form
        linked below to find valid names or ID numbers.
    id_type : `str`
        If 'majorbody', search by name for planets, satellites, or other major bodies.
        If 'smallbody', search by name for asteroids or comets.
        If 'id', search by ID number.
    time : {parse_time_types}
        Time to use in a parse_time-compatible format

    Keyword Arguments
    -----------------
    include_velocity : `bool`
        If True, include the body's velocity in the output coordinate.  Defaults to False.

    Returns
    -------
    `~astropy.coordinates.SkyCoord`
        Location of the solar-system body

    Notes
    -----
    Be aware that there can be discrepancies between the coordinates returned by JPL HORIZONS,
    the coordinates reported in mission data files, and the coordinates returned by
    `~sunpy.coordinates.get_body_heliographic_stonyhurst`.

    References
    ----------
    * `JPL HORIZONS <https://ssd.jpl.nasa.gov/?horizons>`_
    * `JPL HORIZONS form to search bodies <https://ssd.jpl.nasa.gov/horizons.cgi?s_target=1#top>`_
    * `Astroquery <https://astroquery.readthedocs.io/en/latest/>`_

    Examples
    --------
    >>> from sunpy.coordinates.ephemeris import get_horizons_coord

    Query the location of Venus

    >>> get_horizons_coord('Venus barycenter', '2001-02-03 04:05:06')  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for Venus Barycenter (2) [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=2001-02-03T04:05:06.000): (lon, lat, radius) in (deg, deg, AU)
        (-33.93155836, -1.64998443, 0.71915147)>

    Query the location of the SDO spacecraft

    >>> get_horizons_coord('SDO', '2011-11-11 11:11:11')  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for Solar Dynamics Observatory (spac [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=2011-11-11T11:11:11.000): (lon, lat, radius) in (deg, deg, AU)
        (0.01019118, 3.29640728, 0.99011042)>

    Query the location of the SOHO spacecraft via its ID number (-21)

    >>> get_horizons_coord(-21, '2004-05-06 11:22:33', 'id')  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for SOHO (spacecraft) (-21) [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=2004-05-06T11:22:33.000): (lon, lat, radius) in (deg, deg, AU)
        (0.25234902, -3.55863633, 0.99923086)>

    Query the location and velocity of the asteroid Juno

    >>> get_horizons_coord('Juno', '1995-07-18 07:17', 'smallbody', include_velocity=True)  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for 3 Juno (A804 RA) [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=1995-07-18T07:17:00.000): (lon, lat, radius) in (deg, deg, AU)
        (-25.16107532, 14.59098438, 3.17667664)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (-0.03306548, 0.00052415, -2.66709222)>
    """
    obstime = parse_time(time)
    array_time = np.reshape(obstime, (-1,))  # Convert to an array, even if scalar

    # Import here so that astroquery is not a module-level dependency
    from astroquery.jplhorizons import Horizons
    query = Horizons(id=body, id_type=id_type,
                     location='500@10',      # Heliocentric (mean ecliptic)
                     epochs=array_time.tdb.jd.tolist())  # Time must be provided in JD TDB
    try:
        result = query.vectors()
    except Exception as e:  # Catch and re-raise all exceptions, and also provide query URL if generated
        if query.uri is not None:
            log.error(f"See the raw output from the JPL HORIZONS query at {query.uri}")
        raise e
    finally:
        query._session.close()
    log.info(f"Obtained JPL HORIZONS location for {result[0]['targetname']}")
    log.debug(f"See the raw output from the JPL HORIZONS query at {query.uri}")

    # JPL HORIZONS results are sorted by observation time, so this sorting needs to be undone.
    # Calling argsort() on an array returns the sequence of indices of the unsorted list to put the
    # list in order.  Calling argsort() again on the output of argsort() reverses the mapping:
    # the output is the sequence of indices of the sorted list to put that list back in the
    # original unsorted order.
    unsorted_indices = obstime.argsort().argsort()
    result = result[unsorted_indices]

    vector = CartesianRepresentation(result['x'], result['y'], result['z'])
    if include_velocity:
        velocity = CartesianDifferential(result['vx'], result['vy'], result['vz'])
        vector = vector.with_differentials(velocity)
    coord = SkyCoord(vector, frame=HeliocentricEclipticIAU76, obstime=obstime)

    return coord.transform_to(HeliographicStonyhurst).reshape(obstime.shape)
