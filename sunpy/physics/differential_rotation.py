from copy import deepcopy

import numpy as np

import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame, SkyCoord
from astropy.time import TimeDelta

import sunpy.sun.models
from sunpy.coordinates import (
    Heliocentric,
    HeliographicStonyhurst,
    Helioprojective,
    get_earth,
    transform_with_sun_center,
)
from sunpy.coordinates.utils import coordinate_is_on_solar_disk
from sunpy.map import (
    contains_full_disk,
    is_all_off_disk,
    is_all_on_disk,
    map_edges,
    on_disk_bounding_coordinates,
)
from sunpy.map.header_helper import get_observer_meta
from sunpy.time import parse_time
from sunpy.util import expand_list
from sunpy.util.exceptions import warn_user

__all__ = ['solar_rotate_coordinate', 'differential_rotate']


def _validate_observer_args(initial_obstime, observer, time):
    if (observer is not None) and (time is not None):
        raise ValueError(
            "Either the 'observer' or the 'time' keyword must be specified, "
            "but not both simultaneously.")
    elif observer is not None:
        # Check that the new_observer is specified correctly.
        if not (isinstance(observer, BaseCoordinateFrame | SkyCoord)):
            raise ValueError(
                "The 'observer' must be an astropy.coordinates.BaseCoordinateFrame or an astropy.coordinates.SkyCoord.")
        if observer.obstime is None:
            raise ValueError("The observer 'obstime' property must not be None.")
    elif observer is None and time is None:
        raise ValueError("Either the 'observer' or the 'time' keyword must not be None.")


def _get_new_observer(initial_obstime, observer, time):
    """
    Helper function that interprets the possible ways of specifying the
    input to the solar coordinate rotation function.

    If the "observer" argument is not `None`, it is used to specify the location
    of the new observer in space and time.

    If the "time" argument is not `None`, it is used to calculate the duration
    over which to the amount of solar rotation is calculated. Note that using
    the "time" keyword assumes that the new observer is on the Earth. This may
    be a reasonable assumption depending on the application.
    """
    if observer is None:
        obstime = parse_time(time) + initial_obstime
        observer = get_earth(obstime)
    else:
        obstime = observer.obstime
    return observer, obstime


@u.quantity_input
def solar_rotate_coordinate(coordinate, observer=None, time=None, **diff_rot_kwargs):
    """
    Given a coordinate on the Sun, calculate where that coordinate maps to
    at some later time, as a result of solar differential rotation.

    The returned coordinate will have an observation time and observer location
    according to the inputs described below.

    .. note::
        The functions assumes the input coordinate is a position on the Sun
        and therefore ignores the radial component of the 3D input coordinate.

    Parameters
    ----------
    coordinate : `~astropy.coordinates.SkyCoord`
        Any valid coordinate which is transformable to a Stonyhurst Heliographic
        coordinate.

    observer : `~astropy.coordinates.SkyCoord` or None, optional
        The location and observation time for the final observation of the
        rotated coordinate. If the 'time' keyword is specified, this keyword
        must be `None`. The observation time of the output coordinate is
        ``observer.obstime``. When transforming the input coordinate to the
        output coordinate, the intermediate rotation calculations assume the
        observations are made from Earth.

    time : `~sunpy.time.parse_time` compatible object, optional
        The time difference between the start and end time. The keyword
        'time' must be of type `~astropy.time.TimeDelta`,
        `~astropy.units.Quantity` or a value that can be
        parsed using `~astropy.time.TimeDelta`. The keyword is used as the
        duration over which to calculate the apparent motion of a solar feature
        due to differential rotation. The output coordinate will have an
        observation time of the initial coordinate plus this time.

    diff_rot_kwargs : `dict`
        Keyword arguments are passed on to `sunpy.physics.differential_rotation.diff_rot`.

    Returns
    -------
    `~astropy.coordinates.SkyCoord`
        The locations of the input coordinate after differential rotation and
        with observation time and observer according to `time` or `observer`.

    Notes
    -----
    The longitude and latitude columns of the input coordinate influence the
    differential rotation calculation via the definition of the
    'rot' keyword from the `sunpy.physics.differential_rotation.diff_rot`
    function. In the default case, where `time` is provided and `observer` is
    not, the time difference is used to compute the amount of differential
    rotation, i.e. the input `time` keyword is *not* the output time itself.
    The output observation time is the input observation time plus the provided
    `time` keyword.

    Examples
    --------
    A coordinate on the Sun at a specific time

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> from sunpy.coordinates import frames
    >>> from sunpy.physics.differential_rotation import solar_rotate_coordinate
    >>> start_time = '2010-09-10 12:34:56'
    >>> end_time = '2010-09-12 01:00:00'
    >>> c = SkyCoord(-512*u.arcsec, 165*u.arcsec, obstime=start_time, observer="earth", frame=frames.Helioprojective)
    >>> solar_rotate_coordinate(c, time=end_time)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
    <SkyCoord (Helioprojective: obstime=2010-09-12T01:00:00.000, rsun=695700.0 km, observer=<HeliographicStonyhurst Coordinate for 'earth'>): (Tx, Ty, distance) in (arcsec, arcsec, AU)
        (-298.64620102, 166.60837391, 1.00633102)>
    """
    initial_obstime = coordinate.obstime
    observer, new_obstime = _get_new_observer(initial_obstime, observer, time)

    # Calculate the differential rotation in longitude and latitude.
    new_lon = _diff_rot(coordinate, initial_obstime, new_obstime, **diff_rot_kwargs)

    # The input coordinate is transformed to Heliographic Stonyhurst and back to
    # get the correct latitude of the output coordinate (i.e. no latitude
    # dependence on rotational evolution).
    hgc = coordinate.transform_to(HeliographicStonyhurst(obstime=initial_obstime))

    # Return the coordinate as SkyCoord in the HGS frame so that it is easy for
    # users to transform.
    new_coord = SkyCoord(new_lon, hgc.lat, obstime=new_obstime, observer=observer, frame=HeliographicStonyhurst)

    # Convert from HGS to whatever the frame of the input coordinate was
    return new_coord.transform_to(coordinate.frame.replicate_without_data(obstime=new_obstime, observer=observer))


def _transform_to_observer_lonlat(coordinate, obstime):
    """
    Helper function that converts the input coordinate to the observer-centered
    longitude and latitude.
    """
    if not isinstance(coordinate.frame, HeliographicStonyhurst):
        coordinate = coordinate.transform_to(HeliographicStonyhurst(obstime=obstime))
    return coordinate.lon, coordinate.lat


def _find_coordinate_corner_indices(coords):
    """
    Returns the indices of the four corner coordinates of the input coordinate array.

    Parameters
    ----------
    coords : `~astropy.coordinates.SkyCoord`
        Coordinate array of shape "(N, M, ...)".

    Returns
    -------
    `tuple`
        Indices of the four corner coordinates as four (N,) shape arrays.
    """
    # If the input is a scalar coordinate, we return an empty tuple.
    if not coords.shape:
        return ()
    # Otherwise find the indices for the four corners of the coordinate array.
    corners = np.array([[0, 0],
                        [0, coords.shape[1] - 1],
                        [coords.shape[0] - 1, 0],
                        [coords.shape[0] - 1, coords.shape[1] - 1]])
    return tuple(np.split(corners, 4, axis=0))


def _get_diff_rot_factors(coordinate, time_start, time_end, **kwargs):
    """
    Helper function that calculates the factor required to scale the
    output from a differential rotation calculation to the specific time range
    of the input coordinate.

    The differential rotation is calculated for each corner of the input
    coordinate array. The minimum and maximum values are used to account for
    the entire coordinate array.
    """
    corners = _find_coordinate_corner_indices(coordinate)
    if len(corners) == 0:
        # Scalar coordinate
        lon, lat = _transform_to_observer_lonlat(coordinate, time_start)
        return _diff_rot_value(lon, lat, time_start, time_end, **kwargs)
    else:
        corner_drots = []
        for corner in corners:
            lon, lat = _transform_to_observer_lonlat(coordinate[tuple(corner[0])], time_start)
            corner_drots.append(_diff_rot_value(lon, lat, time_start, time_end, **kwargs))
        corner_drots = u.Quantity(corner_drots)
        return u.Quantity([corner_drots.min(), corner_drots.max()])


def _diff_rot_value(lon, lat, time_start, time_end, **kwargs):
    dt = TimeDelta(time_end - time_start).to(u.s).value
    return sun.sun.models.diff_rot_value(lon, lat, dt, **kwargs)


def _diff_rot(coordinate, time_start, time_end, **kwargs):
    """
    Calculate the differential rotation value for the input coordinate.

    Parameters
    ----------
    coordinate : `~astropy.coordinates.SkyCoord`
        Coordinate on the Sun.
    time_start : `~astropy.time.Time`
        Start time for the differential rotation calculation.
    time_end : `~astropy.time.Time`
        End time for the differential rotation calculation.
    kwargs : `dict`
        Keyword arguments that are passed to the differential rotation functions.
    """
    lon, lat = _transform_to_observer_lonlat(coordinate, time_start)
    return _diff_rot_value(lon, lat, time_start, time_end, **kwargs)


def _is_missing_data_invalid(input_map):
    """
    Returns `True` if the missing data value of the input map is not equal to
    the expected integer value expected by the ``WCS.expmap_rotate`` function.
    """
    return input_map.mask is not None and not np.array_equal(input_map.mask, ~np.isfinite(input_map.data))


@u.quantity_input
def differential_rotate(smap, *, time=None, observer=None, **diff_rot_kwargs):
    """
    Given a map, calculate how the map data is transformed as a result of solar
    differential rotation.

    This function will transform a map so that solar features have been
    differentially rotated to their apparent positions at the
    observation time according to the `time` or `observer` inputs, or equivalently
    it will de-rotate a map from the observation time of the map to the observation
    time according to the `time` or `observer` inputs.

    The rotation calculation is performed in the frame of reference of the
    initial observation, and is a simplified calculation of the solar
    differential rotation (i.e., it does not account for variations in
    rotation rate with depth).

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        The input map, which is differentially rotated.
    time : `~sunpy.time.parse_time` compatible object, optional
        The time difference between the start and end time. The keyword
        'time' must be of type `~astropy.time.TimeDelta`,
        `~astropy.units.Quantity` or a value that can be
        parsed using `~astropy.time.TimeDelta`. The keyword is used as the
        duration over which to calculate the apparent motion of a solar feature
        due to differential rotation. The output map will have an observation
        time of the initial map plus this time.
    observer : `~astropy.coordinates.SkyCoord`, optional
        The location and observation time for the final observation of the
        rotated map. The output map will have the same observation time and
        observer as this input. Defaults to None.
        If the 'time' keyword is specified, this keyword
        must be `None`. The observation time of the output map is
        ``observer.obstime``. When transforming the input map to the
        output map, the intermediate rotation calculations assume the
        observations are made from Earth.
    diff_rot_kwargs : `dict`
        Keyword arguments are passed on to `sunpy.physics.differential_rotation.diff_rot`.

    Returns
    -------
    `~sunpy.map.GenericMap`
        The differentially rotated map with the same units as the input map.

    Notes
    -----
    The output map has the same dimensions and data shape as the input map.
    This function does not interpolate from the original map onto a new grid.
    The differential rotation is applied to the image by shifting the
    original image data to new pixel locations. This will cause pixels to be
    wrapped around from one edge of the image to another, and can cause
    artifacts in the output image.

    .. note::

        This function currently only works with maps that have a
        `~sunpy.coordinates.Helioprojective` frame.

    Examples
    --------
    A simple example showing the use of `differential_rotate` for a Helioprojective map.

    >>> import sunpy.map
    >>> from sunpy.data.sample import AIA_171_ROTATED
    >>> import astropy.units as u
    >>> aia_map = sunpy.map.Map(AIA_171_ROTATED)  # doctest: +REMOTE_DATA
    >>> rotated_map = differential_rotate(aia_map, time=12*u.h)  # doctest: +REMOTE_DATA
    """
    # Need to handle the edge case where the entire map is off disk.
    # This is currently a limitation of the approach. We check this first,
    # before any of the transformations take place, so the user gets a clear
    # error message about this.
    if is_all_off_disk(smap):
        warn_user("The entire map is off the solar disk. This can lead to errors in the")

    # Handle maps that are a submap of the full disk (i.e. a zoom-in map)
    # This is a common case, and one that we can handle efficiently.
    # The key is that we want to be careful about the case where the map
    # contains off-disk pixels.
    if not contains_full_disk(smap) and not is_all_on_disk(smap):
        warn_user("Input contains both on-disk and off-disk pixels. This can lead to artifacts.")

    initial_obstime = smap.date
    new_observer, new_obstime = _get_new_observer(initial_obstime, observer, time)

    # Get the differential rotation scaling factor. This is calculated on the
    # corner pixels, and the min/max are used. This is required because the
    # differential rotation is a function of heliographic latitude (and thus
    # varies across the map).
    all_coords = all_coordinates_from_map(smap)
    on_disk = coordinate_is_on_solar_disk(all_coords)

    # Get the edges of the map
    edges = map_edges(smap)
    edge_pixels = np.concatenate(edges).value
    edge_coords = smap.wcs.pixel_to_world(edge_pixels[:, 0], edge_pixels[:, 1])

    # Determine the rotation for each of the edge pixels
    lon, lat = _transform_to_observer_lonlat(edge_coords, initial_obstime)
    edge_drots = _diff_rot_value(lon, lat, initial_obstime, new_obstime, **diff_rot_kwargs)

    # Determine how much to pad the output by, based on how much the edges move.
    rotation_limit = u.Quantity([edge_drots.min(), edge_drots.max()])
    # Convert the rotation in longitude to a shift in pixels (in the x-direction)
    x_offset_in_pixels = (rotation_limit / smap.scale.axis1)
    # Calculate the output padding, rounding up to the next integer pixel.
    x_pad_width = int(np.ceil(np.max(np.abs(x_offset_in_pixels).to_value(u.pix))))
    y_pad_width = 0

    # The on-disk bounding coordinates are needed to find which coordinates we
    # need to input into the differential rotation functions.
    # This is also used in the final calculation.
    on_disk_bl, on_disk_tr = on_disk_bounding_coordinates(smap)
    on_disk_left_edge_coord = deepcopy(on_disk_bl)
    on_disk_left_edge_coord.Ty = smap.center.Ty
    on_disk_right_edge_coord = deepcopy(on_disk_tr)
    on_disk_right_edge_coord.Ty = smap.center.Ty

    # Now pad the input by the pad width and perform the rotation calculation
    # on this padded array.
    # First, calculate the new reference pixel location after padding.
    ny, nx = smap.data.shape
    padded_data = np.pad(smap.data, ((y_pad_width, y_pad_width), (x_pad_width, x_pad_width)),
                         mode='constant', constant_values=np.nan)
    padded_mask = None
    if smap.mask is not None:
        padded_mask = np.pad(smap.mask, ((y_pad_width, y_pad_width), (x_pad_width, x_pad_width)),
                             mode='constant', constant_values=True)
    # Now create the padded map header and meta
    padded_meta = deepcopy(smap.meta)
    # Update the reference pixel for the padded map
    padded_meta['crpix1'] = float(smap.reference_pixel[0].value) + x_pad_width
    padded_meta['crpix2'] = float(smap.reference_pixel[1].value) + y_pad_width
    padded_map = sunpy.map.Map(padded_data, padded_meta, mask=padded_mask)

    # Work out where the output pixels map to in the input grid.
    # The rotation calculation is done in the input frame.
    # We build a grid of indices covering the output image (padded)
    output_ny, output_nx = padded_map.data.shape
    output_x_indices, output_y_indices = np.meshgrid(np.arange(output_nx), np.arange(output_ny))
    output_indices = np.stack([output_x_indices, output_y_indices], axis=-1).reshape(-1, 2)

    # Convert this output grid to world coordinates (on output map WCS).
    output_coords = padded_map.wcs.pixel_to_world(output_indices[:, 0], output_indices[:, 1])

    # Now figure out what input index each output pixel maps to.
    # We do this by:
    # 1. finding the world coordinate of the output pixel
    # 2. working out how much that world coordinate rotated by (to the initial time)
    # 3. finding what input pixel maps to that reverse-rotated world coordinate
    with transform_with_sun_center():
        output_to_initial_lon_diff = _get_diff_rot_factors(
            output_coords, initial_obstime, new_obstime, **diff_rot_kwargs)
    # We need to work out how much each pixel shifted. We take the min/max
    # rotation and convert to a shift in pixels, and just look up the
    # shifted input location.
    # TODO: do this more precisely per-pixel (this is an approximation that
    # does not work well for large rotation periods).
    output_lon_shift_in_pix = (output_to_initial_lon_diff / padded_map.scale.axis1).to_value(u.pix)
    # We just need the absolute value of the shift range for the lookup
    if np.isscalar(output_lon_shift_in_pix):
        output_lon_shift_in_pix = np.array([output_lon_shift_in_pix, output_lon_shift_in_pix])
    # Shift the output indices by the lon shift to get the input indices we actually sample.
    # The strategy is to try sampling at the position minus half the shift and plus half the shift.
    shift_range = np.max(np.abs(output_lon_shift_in_pix))
    half_shift = int(np.ceil(shift_range / 2))

    # Now build the shifted input indices
    input_indices = deepcopy(output_indices.astype(float))
    # Shift by the negative of the shift range (this moves the lookup backwards in time)
    input_indices[:, 0] -= np.linspace(output_lon_shift_in_pix[0],
                                       output_lon_shift_in_pix[-1],
                                       len(input_indices))

    # Now sample the input at these shifted indices.
    # Clip the indices to the range of the input, then sample with round nearest.
    input_x = np.clip(np.rint(input_indices[:, 0]).astype(int), 0, output_nx - 1)
    input_y = np.clip(np.rint(input_indices[:, 1]).astype(int), 0, output_ny - 1)

    # Take the values from the padded input map and reshape back to the output shape.
    rotated_data = padded_map.data[input_y, input_x].reshape(output_ny, output_nx)

    # Do the same for the mask, if present.
    rotated_mask = padded_mask
    if padded_mask is not None:
        rotated_mask = padded_mask[input_y, input_x].reshape(output_ny, output_nx)

    # Now take the output of the transformation and produce the final map.
    # First create the meta of the output (observer update, date, etc.).
    rotated_meta = deepcopy(padded_meta)
    rotated_meta['date-obs'] = new_obstime.isot
    if hasattr(new_observer, 'frame'):
        observer_meta = get_observer_meta(new_observer.frame)
        rotated_meta.update(observer_meta)

    # Build the rotated map
    rotated_map = sunpy.map.Map(rotated_data, rotated_meta, mask=rotated_mask)

    # Crop the output back to the original dimensions, so the user gets back
    # the same size as the input.
    if x_pad_width > 0 or y_pad_width > 0:
        bl = (x_pad_width, y_pad_width) * u.pix
        tr = (x_pad_width + nx - 1, y_pad_width + ny - 1) * u.pix
        rotated_map = rotated_map.submap(bl, top_right=tr)

    return rotated_map
