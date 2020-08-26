"""
Ephemeris calculations using SunPy coordinate frames
"""
import numpy as np

import astropy.units as u
from astropy.constants import c as speed_of_light
from astropy.coordinates import ICRS, HeliocentricEclipticIAU76, SkyCoord, get_body_barycentric
from astropy.coordinates.representation import CartesianRepresentation

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
def get_body_heliographic_stonyhurst(body, time='now', observer=None):
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

    Returns
    -------
    out : `~sunpy.coordinates.frames.HeliographicStonyhurst`
        Location of the solar-system body in the `~sunpy.coordinates.HeliographicStonyhurst` frame

    Notes
    -----
    There is no correction for aberration due to observer motion.  For a body close to the Sun in
    angular direction relative to the observer, the correction can be negligible because the
    apparent location of the body will shift in tandem with the Sun.
    """
    obstime = parse_time(time)

    if observer is None:
        body_icrs = get_body_barycentric(body, obstime)
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

    body_hgs = ICRS(body_icrs).transform_to(HeliographicStonyhurst(obstime=obstime))

    return body_hgs


@add_common_docstring(**_variables_for_parse_time_docstring())
def get_earth(time='now'):
    """
    Return a `~astropy.coordinates.SkyCoord` for the location of the Earth at a specified time in
    the `~sunpy.coordinates.frames.HeliographicStonyhurst` frame.  The longitude will be 0 by definition.

    Parameters
    ----------
    time : {parse_time_types}
        Time to use in a parse_time-compatible format

    Returns
    -------
    out : `~astropy.coordinates.SkyCoord`
        Location of the Earth in the `~sunpy.coordinates.frames.HeliographicStonyhurst` frame
    """
    earth = get_body_heliographic_stonyhurst('earth', time=time)

    # Explicitly set the longitude to 0
    earth = SkyCoord(0*u.deg, earth.lat, earth.radius, frame=earth)

    return earth


@add_common_docstring(**_variables_for_parse_time_docstring())
def get_horizons_coord(body, time='now', id_type='majorbody'):
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
    >>> from sunpy.coordinates import get_horizons_coord

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
    coord = SkyCoord(vector, frame=HeliocentricEclipticIAU76, obstime=obstime)

    return coord.transform_to(HeliographicStonyhurst).reshape(obstime.shape)
