"""
This module provides routines for applying solar rotation functions to
map sequences.
"""
import copy

import numpy as np
from scipy.ndimage import shift

import astropy.units as u

import sunpy.map
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.util.decorators import deprecated

__author__ = 'J. Ireland'

__all__ = ['calculate_solar_rotate_shift', 'mapsequence_solar_derotate']


@deprecated(since='4.0', alternative='`sunkit_image.coalignment.calculate_solar_rotate_shift`')
def calculate_solar_rotate_shift(mc, layer_index=0, **kwargs):
    """
    Calculate the shift that must be applied to each map contained in a mapsequence
    in order to compensate for solar rotation.

    The center of the map is used to calculate the position of each mapsequence
    layer. Shifts are calculated relative to a specified layer in the mapsequence.
    When using this functionality, it is a good idea to check that the shifts
    that were applied to were reasonable and expected. One way of checking this
    is to animate the original mapsequence, animate the derotated mapsequence, and
    compare the differences you see to the calculated shifts. An example use is
    as follows. If you select data from the SDO cutout service, it is common to
    not use the solar tracking implemented by this service. This is because (at
    time of writing) the solar tracking implemented by that service moves the
    image by single pixels at a time. This is not optimal for many use cases,
    as it introduces artificial jumps in the data. So with solar tracking not
    chosen, the selected area is like a window through which you can see the
    Sun rotating underneath.

    Parameters
    ----------
    mc : `sunpy.map.MapSequence`
        The input mapsequence.
    layer_index : int
        The index layer.  Shifts are calculated relative to the time of
        this layer.
    ``**kwargs``
        These keywords are passed to the function
        `sunpy.physics.differential_rotation.solar_rotate_coordinate`.
    Returns
    -------
    x, y : `~astropy.units.Quantity`, ~astropy.units.Quantity`
        The shifts relative to the index layer that can be applied
        to the input mapsequence in order to compensate for solar rotation.
        The shifts are given in arcseconds as understood in helioprojective
        coordinates systems.
    """
    # Size of the data
    nt = len(mc.maps)

    # Storage for the shifts in arcseconds
    xshift_arcseconds = np.zeros(nt) * u.arcsec
    yshift_arcseconds = np.zeros_like(xshift_arcseconds)

    # Layer that
    rotate_to_this_layer = mc.maps[layer_index]

    # Calculate the rotations and the shifts
    for i, m in enumerate(mc):
        # Skip the reference layer
        if i == layer_index:
            continue

        # Calculate the rotation of the center of the map 'm' at its
        # observation time to the observation time of the reference layer
        # indicated by "layer_index".
        new_coordinate = solar_rotate_coordinate(m.center,
                                                 observer=rotate_to_this_layer.observer_coordinate,
                                                 **kwargs)

        # Calculate the shift in arcseconds
        xshift_arcseconds[i] = new_coordinate.Tx - rotate_to_this_layer.center.Tx
        yshift_arcseconds[i] = new_coordinate.Ty - rotate_to_this_layer.center.Ty

    return {"x": xshift_arcseconds, "y": yshift_arcseconds}


@deprecated(since='4.0', alternative='`sunkit_image.coalignment.mapsequence_coalign_by_rotation`')
def mapsequence_solar_derotate(mc, layer_index=0, clip=True, shift=None, **kwargs):
    """
    Move the layers in a mapsequence according to the input shifts.
    If an input shift is not given, the shifts due to
    solar rotation relative to an index layer is calculated and
    applied.  When using this functionality, it is a good idea to check
    that the shifts that were applied to were reasonable and expected.
    One way of checking this is to animate the original mapsequence, animate
    the derotated mapsequence, and compare the differences you see to the
    calculated shifts.

    Parameters
    ----------
    mc : `sunpy.map.MapSequence`
        A mapsequence of shape (ny, nx, nt), where nt is the number of layers in
        the mapsequence.
    layer_index : int
        Solar derotation shifts of all maps in the mapsequence are assumed
        to be relative to the layer in the mapsequence indexed by layer_index.
    clip : bool
        If True, then clip off x, y edges in the datasequence that are potentially
        affected by edges effects.
    ``**kwargs``
        These keywords are passed to the function
        `sunpy.physics.solar_rotation.calculate_solar_rotate_shift`.

    Returns
    -------
    output : `sunpy.map.MapSequence`
        The results of the shifts applied to the input mapsequence.

    Examples
    --------
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> from sunpy.physics.solar_rotation import mapsequence_solar_derotate
    >>> map1 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
    >>> map2 = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE)  # doctest: +REMOTE_DATA
    >>> mc = sunpy.map.Map([map1, map2], sequence=True)  # doctest: +REMOTE_DATA
    >>> derotated_mc = mapsequence_solar_derotate(mc)  # doctest: +SKIP
    >>> derotated_mc = mapsequence_solar_derotate(mc, layer_index=-1)  # doctest: +SKIP
    >>> derotated_mc = mapsequence_solar_derotate(mc, clip=False)  # doctest: +SKIP
    """

    # Size of the data
    nt = len(mc.maps)

    # Storage for the pixel shifts and the shifts in arcseconds
    xshift_keep = np.zeros(nt) * u.pix
    yshift_keep = np.zeros_like(xshift_keep)

    # If no shifts are passed in, calculate them.  Otherwise,
    # use the shifts passed in.
    if shift is None:
        shift = calculate_solar_rotate_shift(mc, layer_index=layer_index, **kwargs)
    xshift_arcseconds = shift['x']
    yshift_arcseconds = shift['y']

    # Calculate the pixel shifts
    for i, m in enumerate(mc):
        xshift_keep[i] = xshift_arcseconds[i] / m.scale[0]
        yshift_keep[i] = yshift_arcseconds[i] / m.scale[1]

    # Apply the pixel shifts and return the mapsequence
    return _apply_shifts(mc, yshift_keep, xshift_keep, clip=clip)


def _apply_shifts(mc, yshift: u.pix, xshift: u.pix, clip=True, **kwargs):
    """
    Apply a set of pixel shifts to a `~sunpy.map.MapSequence`, and return a new
    `~sunpy.map.MapSequence`.

    Parameters
    ----------
    mc : `sunpy.map.MapSequence`
        A `~sunpy.map.MapSequence` of shape ``(ny, nx, nt)``, where ``nt`` is the number of
        layers in the `~sunpy.map.MapSequence`. ``ny`` is the number of pixels in the
        "y" direction, ``nx`` is the number of pixels in the "x" direction.
    yshift : `~astropy.units.Quantity`
        An array of pixel shifts in the y-direction for an image.
    xshift : `~astropy.units.Quantity`
        An array of pixel shifts in the x-direction for an image.
    clip : `bool`, optional
        If `True` (default), then clip off "x", "y" edges of the maps in the sequence that are
        potentially affected by edges effects.

    Notes
    -----
    All other keywords are passed to `scipy.ndimage.shift`.

    Returns
    -------
    `sunpy.map.MapSequence`
        A `~sunpy.map.MapSequence` of the same shape as the input. All layers in
        the `~sunpy.map.MapSequence` have been shifted according the input shifts.
    """
    # New mapsequence will be constructed from this list
    new_mc = []

    # Calculate the clipping
    if clip:
        yclips, xclips = _calculate_clipping(-yshift, -xshift)

    # Shift the data and construct the mapsequence
    for i, m in enumerate(mc):
        shifted_data = shift(copy.deepcopy(m.data), [yshift[i].value, xshift[i].value], **kwargs)
        new_meta = copy.deepcopy(m.meta)
        # Clip if required.  Use the submap function to return the appropriate
        # portion of the data.
        if clip:
            shifted_data = _clip_edges(shifted_data, yclips, xclips)
            new_meta['naxis1'] = shifted_data.shape[1]
            new_meta['naxis2'] = shifted_data.shape[0]
            # Add one to go from zero-based to one-based indexing
            new_meta['crpix1'] = m.reference_pixel.x.value + 1 + xshift[i].value - xshift[0].value
            new_meta['crpix2'] = m.reference_pixel.y.value + 1 + yshift[i].value - yshift[0].value

        new_map = sunpy.map.Map(shifted_data, new_meta)

        # Append to the list
        new_mc.append(new_map)

    return sunpy.map.Map(new_mc, sequence=True)


def _clip_edges(data, yclips: u.pix, xclips: u.pix):
    """
    Clips off the "y" and "x" edges of a 2D array according to a list of pixel
    values. This function is useful for removing data at the edge of 2d images
    that may be affected by shifts from solar de- rotation and layer co-
    registration, leaving an image unaffected by edge effects.

    Parameters
    ----------
    data : `numpy.ndarray`
        A numpy array of shape ``(ny, nx)``.
    yclips : `astropy.units.Quantity`
        The amount to clip in the y-direction of the data. Has units of
        pixels, and values should be whole non-negative numbers.
    xclips : `astropy.units.Quantity`
        The amount to clip in the x-direction of the data. Has units of
        pixels, and values should be whole non-negative numbers.

    Returns
    -------
    `numpy.ndarray`
        A 2D image with edges clipped off according to ``yclips`` and ``xclips``
        arrays.
    """
    ny = data.shape[0]
    nx = data.shape[1]
    # The purpose of the int below is to ensure integer type since by default
    # astropy quantities are converted to floats.
    return data[int(yclips[0].value): ny - int(yclips[1].value),
                int(xclips[0].value): nx - int(xclips[1].value)]


def _calculate_clipping(y: u.pix, x: u.pix):
    """
    Return the upper and lower clipping values for the "y" and "x" directions.

    Parameters
    ----------
    y : `astropy.units.Quantity`
        An array of pixel shifts in the y-direction for an image.
    x : `astropy.units.Quantity`
        An array of pixel shifts in the x-direction for an image.

    Returns
    -------
    `tuple`
        The tuple is of the form ``([y0, y1], [x0, x1])``.
        The number of (integer) pixels that need to be clipped off at each
        edge in an image. The first element in the tuple is a list that gives
        the number of pixels to clip in the y-direction. The first element in
        that list is the number of rows to clip at the lower edge of the image
        in y. The clipped image has "clipping[0][0]" rows removed from its
        lower edge when compared to the original image. The second element in
        that list is the number of rows to clip at the upper edge of the image
        in y. The clipped image has "clipping[0][1]" rows removed from its
        upper edge when compared to the original image. The second element in
        the "clipping" tuple applies similarly to the x-direction (image
        columns). The parameters ``y0, y1, x0, x1`` have the type
        `~astropy.units.Quantity`.
    """
    return ([_lower_clip(y.value), _upper_clip(y.value)] * u.pix,
            [_lower_clip(x.value), _upper_clip(x.value)] * u.pix)


def _upper_clip(z):
    """
    Find smallest integer bigger than all the positive entries in the input
    array.
    """
    zupper = 0
    zcond = z >= 0
    if np.any(zcond):
        zupper = int(np.max(np.ceil(z[zcond])))
    return zupper


def _lower_clip(z):
    """
    Find smallest positive integer bigger than the absolute values of the
    negative entries in the input array.
    """
    zlower = 0
    zcond = z <= 0
    if np.any(zcond):
        zlower = int(np.max(np.ceil(-z[zcond])))
    return zlower
