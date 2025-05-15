"""
Ephemeris calculations using SunPy coordinate frames
"""
import re

import numpy as np
import requests

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
from astropy.io import ascii
from astropy.time import Time

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
def get_body_heliographic_stonyhurst(body, time='now', observer=None, *, include_velocity=False,
                                     quiet=False):
    """
    Return a `~sunpy.coordinates.frames.HeliographicStonyhurst` frame for the location of a
    solar-system body at a specified time. The location can be corrected for light travel time
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
    include_velocity : `bool`, optional
        If True, include the body's velocity in the output coordinate. Defaults to False.
    quiet : `bool`, optional
        If True, the function will not emit logger output for light-travel-time corrections.
        Defaults to False.

    Returns
    -------
    out : `~sunpy.coordinates.frames.HeliographicStonyhurst`
        Location of the solar-system body in the `~sunpy.coordinates.HeliographicStonyhurst` frame

    Notes
    -----
    There is no correction for aberration due to observer motion. For a body close to the Sun in
    angular direction relative to the observer, the correction can be negligible because the
    apparent location of the body will shift in tandem with the Sun.

    For planets other than Earth, Astropy's built-in ephemeris is not as accurate as JPL
    ephemerides, so one can use `astropy.coordinates.solar_system_ephemeris` to switch
    to a JPL ephemeris. See :ref:`astropy-coordinates-solarsystem` for more information, and see
    :ref:`sphx_glr_generated_gallery_units_and_coordinates_venus_transit.py` for an example.

    Examples
    --------
    >>> from sunpy.coordinates.ephemeris import get_body_heliographic_stonyhurst

    Obtain the location of Venus

    >>> get_body_heliographic_stonyhurst('venus', '2012-06-06 04:07:29')
    <HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (0.07349535, 0.05223575, 0.72605496)>

    Obtain the location of Venus as seen from Earth when adjusted for light travel time

    >>> earth = get_body_heliographic_stonyhurst('earth', '2012-06-06 04:07:29')
    >>> get_body_heliographic_stonyhurst('venus', '2012-06-06 04:07:29', observer=earth)
    INFO: Apparent body location accounts for 144.07 seconds of light travel time [sunpy.coordinates.ephemeris]
    <HeliographicStonyhurst Coordinate (obstime=2012-06-06T04:07:29.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (0.07084926, 0.0520573, 0.72605477)>

    Obtain the location and velocity of Mars

    >>> mars = get_body_heliographic_stonyhurst('mars', '2001-02-03', include_velocity=True)
    >>> mars
    <HeliographicStonyhurst Coordinate (obstime=2001-02-03T00:00:00.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (63.03105777, -5.20656151, 1.6251161)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (-0.02323686, 0.00073376, -1.4798387)>
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

        if not quiet:
            if light_travel_time.isscalar:
                ltt_string = f"{light_travel_time.to_value('s'):.2f}"
            else:
                ltt_string = f"{light_travel_time.to_value('s')}"
            log.info(f"Apparent body location accounts for {ltt_string} "
                     "seconds of light travel time")

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
    the `~sunpy.coordinates.frames.HeliographicStonyhurst` frame. The longitude will be zero by
    definition.

    Parameters
    ----------
    time : {parse_time_types}
        Time to use in a parse_time-compatible format
    include_velocity : `bool`, optional
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
    <SkyCoord (HeliographicStonyhurst: obstime=2001-02-03T04:05:06.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (0., -6.18656962, 0.98567647)>
    >>> get_earth('2001-02-03 04:05:06', include_velocity=True)
    <SkyCoord (HeliographicStonyhurst: obstime=2001-02-03T04:05:06.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
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
def get_horizons_coord(body, time='now', id_type=None, *, include_velocity=False):
    """
    Queries JPL HORIZONS and returns a `~astropy.coordinates.SkyCoord` for the location of a
    solar-system body at a specified time. This location is the instantaneous or "true" location,
    and is not corrected for light travel time or observer motion.

    Parameters
    ----------
    body : `str`
        The solar-system body for which to calculate positions. One can also use the search form
        linked below to find valid names or ID numbers.
    id_type : `None`, `str`
        Defaults to `None`, which searches major bodies first, and then searches
        small bodies (comets and asteroids) if no major body is found. If
        ``'smallbody'``, the search is limited to only small bodies. If
        ``'designation'``, the search is limited to only small-body designations.
        If ``'name'``, the search is limited to only small-body names. If
        ``'asteroid_name'`` or ``'comet_name'``, the search is limited to only
        asteroid names or only comet names, respectively.
    time : {parse_time_types}, `dict`
        Time to use in a parse_time-compatible format.

        Alternatively, this can be a dictionary defining a range of times and
        dates; the range dictionary has to be of the form
        {{'start': start_time, 'stop': stop_time, 'step':'n[y|d|m]'}}.
        ``start_time`` and ``stop_time`` must be in a parse_time-compatible format,
        and are interpreted as UTC time. ``step`` must be a string with either a
        number and interval length (e.g. for every 10 minutes, ``'10m'``), or a
        plain number for a number of evenly spaced intervals.

    include_velocity : `bool`, optional
        If True, include the body's velocity in the output coordinate. Defaults to False.

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
    * `JPL HORIZONS <https://ssd.jpl.nasa.gov/?horizons>`__
    * `JPL HORIZONS form to search bodies <https://ssd.jpl.nasa.gov/horizons.cgi?s_target=1#top>`__

    Examples
    --------
    .. minigallery:: sunpy.coordinates.get_horizons_coord

    >>> from sunpy.coordinates import get_horizons_coord

    Query the location of Venus

    >>> get_horizons_coord('Venus barycenter', '2001-02-03 04:05:06')  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for Venus Barycenter (2) [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=2001-02-03T04:05:06.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (-33.93155836, -1.64998443, 0.71915147)>

    Query the location of the SDO spacecraft

    >>> get_horizons_coord('SDO', '2011-11-11 11:11:11')  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for Solar Dynamics Observatory (spacecraft) (-136395) [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=2011-11-11T11:11:11.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (0.01019118, 3.29640728, 0.99011042)>

    Query the location of the SOHO spacecraft via its ID number (-21)

    >>> get_horizons_coord(-21, '2004-05-06 11:22:33')  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for SOHO (spacecraft) (-21) [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=2004-05-06T11:22:33.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (0.25234902, -3.55863633, 0.99923086)>

    Query the location and velocity of the asteroid Juno

    >>> get_horizons_coord('Juno', '1995-07-18 07:17', 'smallbody', include_velocity=True)  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for 3 Juno (A804 RA) [sunpy.coordinates.ephemeris]
    <SkyCoord (HeliographicStonyhurst: obstime=1995-07-18T07:17:00.000, rsun=695700.0 km): (lon, lat, radius) in (deg, deg, AU)
        (-25.16107532, 14.59098438, 3.17667664)
     (d_lon, d_lat, d_radius) in (arcsec / s, arcsec / s, km / s)
        (-0.03306548, 0.00052415, -2.66709222)>

    Query the location of Solar Orbiter at a set of 12 regularly sampled times

    >>> get_horizons_coord('Solar Orbiter',
    ...                    time={{'start': '2020-12-01',
    ...                           'stop': '2020-12-02',
    ...                           'step': '12'}})  # doctest: +REMOTE_DATA
    INFO: Obtained JPL HORIZONS location for Solar Orbiter (spacecraft) (-144) [sunpy.coordinates.ephemeris]
    ...
    """
    # Reference plane defaults to the ecliptic (IAU 1976 definition)
    args = {
        'EPHEM_TYPE': 'VECTORS',
        'OUT_UNITS': 'AU-D',  # units of AU and days
        'CENTER': '500@10',  # origin is body center (500) of the Sun (10)
        'VEC_TABLE': '2',  # output the 6-element state vector
        'CSV_FORMAT': 'YES',
    }

    if id_type in [None, 'smallbody', 'designation', 'name', 'asteroid_name', 'comet_name']:
        prepend = {
            None: "",
            'smallbody': "",
            'designation': "DES=",
            'name': "NAME=",
            'asteroid_name': "ASTNAM=",
            'comet_name': "COMNAM=",
        }
        if id_type == 'smallbody' and body[-1] != ';':
            body += ';'
        args['COMMAND'] = f"'{prepend[id_type]}{body}'"
    else:
        raise ValueError("Invalid id_type")

    if isinstance(time, dict):
        if set(time.keys()) != set(['start', 'stop', 'step']):
            raise ValueError('time dictionary must have the keys ["start", "stop", "step"]')
        jpl_fmt = "'%Y-%m-%d %H:%M:%S.%f'"
        args['START_TIME'] = parse_time(time['start']).tdb.strftime(jpl_fmt)
        args['STOP_TIME'] = parse_time(time['stop']).tdb.strftime(jpl_fmt)
        args['STEP_SIZE'] = time['step']
    else:
        obstime = parse_time(time)
        if obstime.size >= 10000:
            raise ValueError("For more than 10,000 time values, use dictionary input.")
        array_time = np.reshape(obstime, (-1,))  # Convert to an array, even if scalar
        epochs = array_time.tdb.jd.tolist()  # Time must be provided in JD TDB
        args['TLIST'] = '\n'.join(str(epoch) for epoch in epochs)
        args['TLIST_TYPE'] = 'JD'  # needs to be explicitly set to avoid potential MJD confusion

    contents = "!$$SOF\n" + '\n'.join(f"{k}={v}" for k, v in args.items())
    log.debug(f"JPL HORIZONS query via POST request:\n{contents}")

    output = requests.post('https://ssd.jpl.nasa.gov/api/horizons_file.api',
                           data={'format': 'text'}, files={'input': contents})

    lines = output.text.splitlines()
    start_index, stop_index = 0, len(lines)
    success = False
    error_message = ''
    for index, line in enumerate(lines):
        if line.startswith("Target body name:"):
            target_name = re.search(r': (.*) {', line).group(1)
            log.info(f"Obtained JPL HORIZONS location for {target_name}")

        if "Multiple major-bodies match string" in line:
            error_message = '\n'.join(lines[index:-1])
            break

        if "Matching small-bodies:" in line:
            # Prepare the error message assuming there are multiple matches
            error_message = '\n'.join(lines[index:-1])

            # Change the error message if there are actually zero matches
            if "No matches found." in error_message:
                error_message = "No matches found."

            break

        if line.startswith("$$SOE"):
            start_index = index + 1
        if line.startswith("$$EOE"):
            stop_index = index
            success = True

    if not success:
        if error_message:
            raise ValueError(error_message)
        else:
            raise RuntimeError(f"Unknown JPL HORIZONS error:\n{output.text}")

    column_names = [name.strip() for name in lines[start_index - 3].split(',')]
    result = ascii.read(lines[start_index:stop_index], names=column_names)

    if isinstance(time, dict):
        obstime_tdb = parse_time(result['JDTDB'], format='jd', scale='tdb')
        obstime = Time(obstime_tdb, format='isot', scale='utc')
    else:
        # JPL HORIZONS results are sorted by observation time, so this sorting needs to be undone.
        # Calling argsort() on an array returns the sequence of indices of the unsorted list to put the
        # list in order. Calling argsort() again on the output of argsort() reverses the mapping:
        # the output is the sequence of indices of the sorted list to put that list back in the
        # original unsorted order.
        unsorted_indices = obstime.argsort().argsort()
        result = result[unsorted_indices]

    vector = CartesianRepresentation(result['X'], result['Y'], result['Z']) * u.AU
    if include_velocity:
        velocity = CartesianDifferential(result['VX'], result['VY'], result['VZ']) * u.AU / u.day
        vector = vector.with_differentials(velocity)
    coord = SkyCoord(vector, frame=HeliocentricEclipticIAU76, obstime=obstime)

    return coord.transform_to(HeliographicStonyhurst).reshape(obstime.shape)
