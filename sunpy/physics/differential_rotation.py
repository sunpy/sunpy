import astropy.units as u

import sunpy.sun.models
from sunpy.coordinates.utils import solar_coordinate_rotation
from sunpy.map.maputils import diff_rotation
from sunpy.util.decorators import deprecated

__all__ = ['diff_rot', 'solar_rotate_coordinate', 'differential_rotate']


@u.quantity_input
@deprecated(since="6.0", alternative="sunpy.sun.models.differential_rotation")
def diff_rot(duration: u.s, latitude: u.deg, rot_type='howard', frame_time='sidereal'):
    r"""
    This function computes the change in longitude over days in degrees.

    Parameters
    ----------
    duration : `~astropy.units.Quantity`
        Number of seconds to rotate over.
    latitude : `~astropy.units.Quantity`
        heliographic coordinate latitude in Degrees.
    rot_type : `str`
        The differential rotation model to use.

        One of:

        | ``howard`` : Use values from :cite:t:`howard_solar_1990`
        | ``snodgrass`` : Use values from :cite:t:`snodgrass_magnetic_1983`
        | ``allen`` : Use values from Allen's Astrophysical Quantities, and simpler equation.
        | ``rigid`` : Use values from `~sunpy.sun.constants.sidereal_rotation_rate`.

    frame_time : `str`
        One of : ``'sidereal'`` or  ``'synodic'``. Choose 'type of day' time reference frame.

    Returns
    -------
    longitude_delta : `~astropy.units.Quantity`
        The change in longitude over days (units=degrees)

    Notes
    -----
    The rotation rate at a heliographic latitude :math:`\theta` is given by

    .. math::

        A + B \sin^{2} \left (\theta \right ) + C \sin^{4} \left ( \theta \right )

    where :math:`A, B, C` are constants that depend on the model:

    ========= ======= ====== ====== ==========
    Model     A       B      C      Unit
    ========= ======= ====== ====== ==========
    howard    2.894   -0.428 -0.370 microrad/s
    snodgrass 2.851   -0.343 -0.474 microrad/s
    allen     14.44   -3.0   0      deg/day
    rigid     14.1844 0      0      deg/day
    ========= ======= ====== ====== ==========

    1 microrad/s is approximately 4.95 deg/day.
    See also the comparisons in :cite:t:`beck_comparison_2000`.
    """

    return sunpy.sun.models.differential_rotation(duration, latitude, model=rot_type, frame_time=frame_time)


@deprecated(since="6.1", alternative="sunpy.coordinates.utils.solar_coordinate_rotation")
def solar_rotate_coordinate(coordinate, observer=None, time=None, **diff_rot_kwargs):
    """
    Given a coordinate on the Sun, calculate where that coordinate maps to
    as seen by a new observer at some later or earlier time, given that
    the input coordinate rotates according to the solar rotation profile.

    The amount of solar rotation is based on the amount of time between the
    observation time of the input coordinate and the observation time of the
    new observer. The new observer is specified in one of two ways, either
    using the "observer" or "time" keywords.

    If the "observer" keyword is set, it is used to specify the location
    of the new observer in space and time. The difference between the
    coordinate time and the new observer time is used to calculate the amount
    of solar rotation applied, and the location of the new observer in space
    is used to calculate where the rotated coordinate is as seen from the
    new observer.

    If the "time" keyword is set, it is used to specify the number of
    seconds to rotate the coordinate by. Note that using the "time" keyword
    assumes that the new observer is on the Earth. This may be a reasonable
    assumption depending on the application.

    Either the "observer" or "time" keyword must be specified, but both
    cannot be specified at the same time.

    Parameters
    ----------
    coordinate : `~astropy.coordinates.SkyCoord`
        Any valid coordinate which is transformable to Heliographic Stonyhurst.
    observer : `~astropy.coordinates.BaseCoordinateFrame`, `~astropy.coordinates.SkyCoord`, None
        The location of the new observer in space and time (the observer must have an
        interpretable obstime property).
    time : `~astropy.time.Time`, `~astropy.time.TimeDelta`, `~astropy.units.Quantity`, None
    **diff_rot_kwargs : `dict`
        Keyword arguments are passed on as keyword arguments to `~sunpy.sun.models.differential_rotation`.
        Note that the keyword "frame_time" is automatically set to the value
        "sidereal".

    Returns
    -------
    coordinate : `~astropy.coordinates.SkyCoord`
        The locations of the input coordinates after the application of
        solar rotation as seen from the point-of-view of the new observer.

    Notes
    -----
    The translational motion of the Sun over the time interval will be ignored.
    See :func:`~sunpy.coordinates.transform_with_sun_center`.
    """
    return solar_coordinate_rotation(coordinate, observer, time, **diff_rot_kwargs)


@deprecated(since="6.1", alternative="sunpy.map.maputils.diff_rotation")
def differential_rotate(smap, observer=None, time=None, **diff_rot_kwargs):
    """
    Warp a `~sunpy.map.GenericMap` to take into account both solar differential
    rotation and the changing location of the observer.

    .. warning::
        This function, while greatly improved in 1.0, is still experimental.
        Please validate that it gives you results you expect and report any
        discrepancies on the SunPy issue tracker.

    The function transforms the input map data pixels by first rotating each
    pixel according to solar differential rotation. The amount of solar
    differential applied is calculated by the time difference between the
    observation time of map and the new observation time, as specified by either the
    "time" keyword or the "obstime" property of the "observer" keyword.
    The location of the rotated pixels are then transformed to locations on the Sun
    as seen from the new observer position. This is desirable since in most cases
    the observer does not remain at a fixed position in space. If
    the "time" keyword is used then the new observer position is assumed to
    be based on the location of the Earth. If the "observer" keyword is used then
    this defines the new observer position.

    The function works with full disk maps and maps that contain portions of the
    solar disk (maps that are entirely off-disk will raise an error). When the
    input map contains the full disk, the output map has the same dimensions as
    the input map. When the input map images only part of the solar disk, only
    the on-disk pixels are differentially rotated and the output map can have
    a different dimensions compared to the input map. In this case any off-disk
    emission shown in the input map is not included in the output map.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        Original map that we want to transform.
    observer : `~astropy.coordinates.BaseCoordinateFrame`, `~astropy.coordinates.SkyCoord`, `None`, optional
        The location of the new observer.
        Instruments in Earth orbit can be approximated by using the position
        of the Earth at the observation time of the new observer.
    time : sunpy-compatible time, `~astropy.time.TimeDelta`, `~astropy.units.Quantity`, `None`, optional
        Used to define the duration over which the amount of solar rotation is
        calculated. If 'time' is an `~astropy.time.Time` then the time interval
        is difference between 'time' and the map observation time. If 'time' is
        `~astropy.time.TimeDelta` or `~astropy.units.Quantity` then the calculation
        is "initial_obstime + time".

    Returns
    -------
    `~sunpy.map.GenericMap`
        A map with the result of applying solar differential rotation to the
        input map.

    Notes
    -----
    The translational motion of the Sun over the time interval will be ignored.
    See :func:`~sunpy.coordinates.transform_with_sun_center`.
    """
    return diff_rotation(smap, observer, time, **diff_rot_kwargs)
