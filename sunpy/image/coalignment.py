"""
This module provides routines for the coalignment of images and mapsequences.

Currently this module provides image coalignment by template matching.
Which is partially inspired by the SSWIDL routine
`tr_get_disp.pro <http://www.heliodocs.com/php/xdoc_print.php?file=$SSW/trace/idl/util/tr_get_disp.pro>`_.

In this implementation, the template matching is handled via the scikit-image
routine :func:`skimage.feature.match_template`.

References
----------
Template matching algorithm:

 * http://scribblethink.org/Work/nvisionInterface/nip.html
 * J.P. Lewis, Fast Template Matching, Vision Interface 95, Canadian Image
   Processing and Pattern Recognition Society, Quebec City, Canada, May 15-19,
   1995, p. 120-123 http://www.scribblethink.org/Work/nvisionInterface/vi95_lewis.pdf.
"""
from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.ndimage.interpolation import shift
from copy import deepcopy
from astropy import units as u
# Image co-registration by matching templates
from skimage.feature import match_template

# SunPy imports
import sunpy.map
from sunpy.util import deprecated
from sunpy.map.mapbase import GenericMap

__author__ = 'J. Ireland'

__all__ = ['calculate_shift', 'clip_edges', 'calculate_clipping',
           'match_template_to_layer', 'find_best_match_location',
           'get_correlation_shifts', 'parabolic_turning_point',
           'repair_image_nonfinite', 'apply_shifts',
           'mapcube_coalign_by_match_template',
           'mapsequence_coalign_by_match_template',
           'calculate_match_template_shift']


def _default_fmap_function(data):
    """
    This function ensures that the data are floats.  It is the default data
    manipulation function for the coalignment method.
    """
    return np.float64(data)


def calculate_shift(this_layer, template):
    """Calculates the pixel shift required to put the template in the "best"
    position on a layer.

    Parameters
    ----------
    this_layer : `~numpy.ndarray`
        A numpy array of size (ny, nx), where the first two dimensions are
        spatial dimensions.
    template : `~numpy.ndarray`
        A numpy array of size (N, M) where N < ny and M < nx.

    Returns
    -------
    shifts : tuple
        Pixel shifts (yshift, xshift) relative to the offset of the template
        to the input array.
    """
    # Repair any NANs, Infs, etc in the layer and the template
    this_layer = repair_image_nonfinite(this_layer)
    template = repair_image_nonfinite(template)

    # Calculate the correlation array matching the template to this layer
    corr = match_template_to_layer(this_layer, template)

    # Calculate the y and x shifts in pixels
    return find_best_match_location(corr)


#
# Remove the edges of a datacube
#
@u.quantity_input(yclips=u.pix, xclips=u.pix)
def clip_edges(data, yclips, xclips):
    """
    Clips off the y and x edges of a 2d array according to a list of pixel
    values.  This function is useful for removing data at the edge of
    2d images that may be affected by shifts from solar de-rotation and
    layer co-registration, leaving an image unaffected by edge effects.

    Parameters
    ----------
    data : `numpy.ndarray`
        A numpy array of shape (ny, nx).
    yclips : `astropy.units.Quantity`
        The amount to clip in the y-direction of the data.  Has units of
        pixels, and values should be whole non-negative numbers.
    xclips : `astropy.units.Quantity`
        The amount to clip in the x-direction of the data.  Has units of
        pixels, and values should be whole non-negative numbers.

    Returns
    -------
    image : `numpy.ndarray`
        A 2d image with edges clipped off according to yclips and xclips
        arrays.
    """
    ny = data.shape[0]
    nx = data.shape[1]
    # The purpose of the int below is to ensure integer type since by default
    # astropy quantities are converted to floats.
    return data[int(yclips[0].value): ny - int(yclips[1].value),
                int(xclips[0].value): nx - int(xclips[1].value)]


#
# Return the upper and lower clipping values for the y and x directions an
# input set of pixel shifts y and x
#
@u.quantity_input(y=u.pix, x=u.pix)
def calculate_clipping(y, x):
    """
    Return the upper and lower clipping values for the y and x directions.

    Parameters
    ----------
    y : `astropy.units.Quantity`
        An array of pixel shifts in the y-direction for an image.
    x : `astropy.units.Quantity`
        An array of pixel shifts in the x-direction for an image.

    Returns
    -------
    clipping : tuple
        The tuple is of the form ([y0, y1], [x0, x1]).
        The number of (integer) pixels that need to be clipped off at each
        edge in an image. The first element in the tuple is a list that gives
        the number of pixels to clip in the y-direction.  The first element in
        that list is the number of rows to clip at the lower edge of the image
        in y.  The clipped image has "clipping[0][0]" rows removed from its
        lower edge when compared to the original image.  The second element in
        that list is the number of rows to clip at the upper edge of the image
        in y.  The clipped image has "clipping[0][1]" rows removed from its
        upper edge when compared to the original image.  The second element in
        the "clipping" tuple applies similarly to the x-direction (image
        columns).  The parameters y0, y1, x0, x1 have the type
        `~astropy.units.Quantity`.
    """
    return ([_lower_clip(y.value), _upper_clip(y.value)] * u.pix,
            [_lower_clip(x.value), _upper_clip(x.value)] * u.pix)


#
# Helper functions for clipping edges
#
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


def match_template_to_layer(layer, template):
    """
    Calculate the correlation array that describes how well the template
    matches the layer. All inputs are assumed to be numpy arrays.  This
    function requires the "match_template" function in scikit image.

    Parameters
    ----------
    layer : `~numpy.ndarray`
        A numpy array of size (ny, nx).
    template : `~numpy.ndarray`
        A numpy array of size (N, M) where N < ny and M < nx.

    Returns
    -------
    correlationarray : `~numpy.ndarray`
        A correlation array between the layer and the template.
        The values in the array range between 0 and 1.
    """
    return match_template(layer, template)


def find_best_match_location(corr):
    """
    Calculate an estimate of the location of the peak of the correlation
    result in image pixels.

    Parameters
    ----------
    corr : `~numpy.ndarray`
        A 2-d correlation array.

    Returns
    -------
    shift : `~astropy.units.Quantity`
        The shift amounts (y, x) in image pixels.  Subpixel values are
        possible.
    """
    # Get the index of the maximum in the correlation function
    ij = np.unravel_index(np.argmax(corr), corr.shape)
    cor_max_x, cor_max_y = ij[::-1]

    # Get the correlation function around the maximum
    array_around_maximum = corr[np.max([0, cor_max_y - 1]): np.min([cor_max_y + 2, corr.shape[0] - 1]),
                                  np.max([0, cor_max_x - 1]): np.min([cor_max_x + 2, corr.shape[1] - 1])]
    y_shift_relative_to_maximum, x_shift_relative_to_maximum = \
    get_correlation_shifts(array_around_maximum)

    # Get shift relative to correlation array
    y_shift_relative_to_correlation_array = y_shift_relative_to_maximum + cor_max_y * u.pix
    x_shift_relative_to_correlation_array = x_shift_relative_to_maximum + cor_max_x * u.pix

    return y_shift_relative_to_correlation_array, x_shift_relative_to_correlation_array


def get_correlation_shifts(array):
    """
    Estimate the location of the maximum of a fit to the input array.  The
    estimation in the x and y directions are done separately. The location
    estimates can be used to implement subpixel shifts between two different
    images.

    Parameters
    ----------
    array : `~numpy.ndarray`
        An array with at least one dimension that has three elements.  The
        input array is at most a 3 x 3 array of correlation values calculated
        by matching a template to an image.

    Returns
    -------
    peakloc : `~astropy.units.Quantity`
        The (y, x) location of the peak of a parabolic fit, in image pixels.
    """
    # Check input shape
    ny = array.shape[0]
    nx = array.shape[1]
    if nx > 3 or ny > 3:
        print('Input array is too big in at least one dimension. Returning Nones')
        return None, None

    # Find where the maximum of the input array is
    ij = np.unravel_index(np.argmax(array), array.shape)
    x_max_location, y_max_location = ij[::-1]

    # Estimate the location of the parabolic peak if there is enough data.
    # Otherwise, just return the location of the maximum in a particular
    # direction.
    if ny == 3:
        y_location = parabolic_turning_point(array[:, x_max_location])
    else:
        y_location = 1.0 * y_max_location

    if nx == 3:
        x_location = parabolic_turning_point(array[y_max_location, :])
    else:
        x_location = 1.0 * x_max_location

    return y_location * u.pix, x_location * u.pix


def parabolic_turning_point(y):
    """
    Find the location of the turning point for a parabola
    y(x) = ax^2 + bx + c, given input values y(-1), y(0), y(1).
    The maximum is located at x0 = -b / 2a .  Assumes
    that the input array represents an equally spaced sampling at the
    locations y(-1), y(0) and y(1).

    Parameters
    ----------
    y : `~numpy.ndarray`
        A one dimensional numpy array of shape 3 with entries that sample the
        parabola at -1, 0, and 1.

    Returns
    -------
    location : float
        A digit, the location of the parabola maximum.
    """
    numerator = -0.5 * y.dot([-1, 0, 1])
    denominator = y.dot([1, -2, 1])
    return numerator / denominator


def repair_image_nonfinite(image):
    """
    Return a new image in which all the nonfinite entries of the original
    image have been replaced by the local mean.

    Parameters
    ----------
    image : `~numpy.ndarray`
        A two-dimensional `~numpy.ndarray`.

    Returns
    -------
    repaired_image : `~numpy.ndarray`
        A two-dimensional `~numpy.ndarray` of the same shape as the input
        that has all the non-finite entries replaced by a local mean.  The
        algorithm repairs one non-finite entry at every pass.  At each pass,
        the next non-finite value is replaced by the mean of its finite
        valued nearest neighbours.
    """
    repaired_image = deepcopy(image)
    nx = repaired_image.shape[1]
    ny = repaired_image.shape[0]
    bad_index = np.where(np.logical_not(np.isfinite(repaired_image)))
    while bad_index[0].size != 0:
        by = bad_index[0][0]
        bx = bad_index[1][0]

        # x locations taking in to account the boundary
        x = bx
        if bx == 0:
            x = 1
        if bx == nx - 1:
            x = nx - 2

        # y locations taking in to account the boundary
        y = by
        if by == 0:
            y = 1
        if by == ny - 1:
            y = ny - 2

        # Get the sub array around the bad index, and find the local mean
        # ignoring nans
        subarray = repaired_image[y - 1: y + 2, x - 1: x + 2]
        repaired_image[by, bx] = np.mean(subarray[np.isfinite(subarray)])
        bad_index = np.where(np.logical_not(np.isfinite(repaired_image)))
    return repaired_image


@u.quantity_input(yshift=u.pix, xshift=u.pix)
def apply_shifts(mc, yshift, xshift, clip=True, **kwargs):
    """
    Apply a set of pixel shifts to a `~sunpy.map.MapSequence`, and return a new
    `~sunpy.map.MapSequence`.

    Parameters
    ----------
    mc : `sunpy.map.MapSequence`
        A `~sunpy.map.MapSequence` of shape (ny, nx, nt), where nt is the number of
        layers in the `~sunpy.map.MapSequence`.  'ny' is the number of pixels in the
        y direction, 'nx' is the number of pixels in the 'x' direction.

    yshift : `~astropy.units.Quantity` instance
        An array of pixel shifts in the y-direction for an image.

    xshift : `~astropy.units.Quantity` instance
        An array of pixel shifts in the x-direction for an image.

    clip : bool
        If True, then clip off x, y edges of the maps in the sequence that are
        potentially affected by edges effects.

    All other keywords are passed to `scipy.ndimage.interpolation.shift`.

    Returns
    -------
    newmapsequence : `sunpy.map.MapSequence`
        A `~sunpy.map.MapSequence` of the same shape as the input.  All layers in
        the `~sunpy.map.MapSequence` have been shifted according the input shifts.
    """
    # New mapsequence will be constructed from this list
    new_mc = []

    # Calculate the clipping
    if clip:
        yclips, xclips = calculate_clipping(-yshift, -xshift)

    # Shift the data and construct the mapsequence
    for i, m in enumerate(mc):
        shifted_data = shift(deepcopy(m.data), [yshift[i].value, xshift[i].value], **kwargs)
        new_meta = deepcopy(m.meta)
        # Clip if required.  Use the submap function to return the appropriate
        # portion of the data.
        if clip:
            shifted_data = clip_edges(shifted_data, yclips, xclips)
            new_meta['naxis1'] = shifted_data.shape[1]
            new_meta['naxis2'] = shifted_data.shape[0]
            new_meta['crpix1'] = m.reference_pixel.x.value + xshift[i].value - xshift[0].value
            new_meta['crpix2'] = m.reference_pixel.y.value + yshift[i].value - yshift[0].value

        new_map = sunpy.map.Map(shifted_data, new_meta)

        # Append to the list
        new_mc.append(new_map)

    return sunpy.map.Map(new_mc, sequence=True)


def calculate_match_template_shift(mc, template=None, layer_index=0,
                                   func=_default_fmap_function):
    """
    Calculate the arcsecond shifts necessary to co-register the layers in a
    `~sunpy.map.MapSequence` according to a template taken from that
    `~sunpy.map.MapSequence`.

    When using this functionality, it is a good idea to check that the shifts
    that were applied to were reasonable and expected. One way of checking this
    is to animate the original `~sunpy.map.MapSequence`, animate the coaligned
    `~sunpy.map.MapSequence`, and compare the differences you see to the
    calculated shifts.

    Parameters
    ----------
    mc : `sunpy.map.MapSequence`
        A `~sunpy.map.MapSequence` of shape (ny, nx, nt), where nt is the number of
        layers in the `~sunpy.map.MapSequence`.

    template : {None | `~sunpy.map.Map` | `~numpy.ndarray`}
        The template used in the matching.  If an `~numpy.ndarray` is passed,
        the `~numpy.ndarray` has to have two dimensions.

    layer_index : int
        The template is assumed to refer to the map in the `~sunpy.map.MapSequence`
        indexed by the value of "layer_index".  Displacements of all maps in the
        `~sunpy.map.MapSequence` are assumed to be relative to this layer.  The
        displacements of the template relative to this layer are therefore
        (0, 0).

    func : function
        A function which is applied to the data values before the coalignment
        method is applied.  This can be useful in coalignment, because it is
        sometimes better to co-align on a function of the data rather than the
        data itself.  The calculated shifts are applied to the original data.
        Examples of useful functions to consider for EUV images are the
        logarithm or the square root.  The function is of the form
        func = F(data).  The default function ensures that the data are
        floats.

    """

    # Size of the data
    ny = mc.maps[layer_index].data.shape[0]
    nx = mc.maps[layer_index].data.shape[1]
    nt = len(mc.maps)

    # Calculate a template.  If no template is passed then define one
    # from the index layer.
    if template is None:
        tplate = mc.maps[layer_index].data[int(ny/4): int(3*ny/4),
                                           int(nx/4): int(3*nx/4)]
    elif isinstance(template, GenericMap):
        tplate = template.data
    elif isinstance(template, np.ndarray):
        tplate = template
    else:
        raise ValueError('Invalid template.')

    # Apply the function to the template
    tplate = func(tplate)

    # Storage for the pixel shift
    xshift_keep = np.zeros(nt) * u.pix
    yshift_keep = np.zeros_like(xshift_keep)

    # Storage for the arcsecond shift
    xshift_arcseconds = np.zeros(nt) * u.arcsec
    yshift_arcseconds = np.zeros_like(xshift_arcseconds)

    # Match the template and calculate shifts
    for i, m in enumerate(mc.maps):
        # Get the next 2-d data array
        this_layer = func(m.data)

        # Calculate the y and x shifts in pixels
        yshift, xshift = calculate_shift(this_layer, tplate)

        # Keep shifts in pixels
        yshift_keep[i] = yshift
        xshift_keep[i] = xshift

    # Calculate shifts relative to the template layer
    yshift_keep = yshift_keep - yshift_keep[layer_index]
    xshift_keep = xshift_keep - xshift_keep[layer_index]

    for i, m in enumerate(mc.maps):
        # Calculate the shifts required in physical units, which are
        # presumed to be arcseconds.
        xshift_arcseconds[i] = xshift_keep[i] * m.scale[0]
        yshift_arcseconds[i] = yshift_keep[i] * m.scale[1]

    return {"x": xshift_arcseconds, "y": yshift_arcseconds}

@deprecated('0.9.1', alternative='mapsequence_coalign_by_match_template')
def mapcube_coalign_by_match_template(mc, template=None, layer_index=0,
                                      func=_default_fmap_function, clip=True,
                                      shift=None, **kwargs):
    """
    Co-register the layers in a `~sunpy.map.MapCube` according to a template
    taken from that `~sunpy.map.MapCube`.

    When using this functionality, it is a good idea to check that the shifts
    that were applied to were reasonable and expected. One way of checking this
    is to animate the original `~sunpy.map.MapCube`, animate the coaligned
    `~sunpy.map.MapCube`, and compare the differences you see to the calculated
    shifts.


    Parameters
    ----------
    mc : `sunpy.map.MapCube`
        A `~sunpy.map.MapCube` of shape (ny, nx, nt), where nt is the number of
        layers in the `~sunpy.map.MapCube`.
    template : {None | sunpy.map.Map | `~numpy.ndarray`}
        The template used in the matching.  If an `~numpy.ndarray` is passed,
        the `~numpy.ndarray` has to have two dimensions.
    layer_index : int
        The template is assumed to refer to the map in the `~sunpy.map.MapCube`
        indexed by the value of "layer_index".  Displacements of all maps in the
        `~sunpy.map.MapCube` are assumed to be relative to this layer.  The
        displacements of the template relative to this layer are therefore
        (0, 0).
    func : function
        A function which is applied to the data values before the coalignment
        method is applied.  This can be useful in coalignment, because it is
        sometimes better to co-align on a function of the data rather than the
        data itself.  The calculated shifts are applied to the original data.
        Examples of useful functions to consider for EUV images are the
        logarithm or the square root.  The function is of the form
        func = F(data).  The default function ensures that the data are
        floats.
    clip : bool
        If True, then clip off x, y edges in the datacube that are potentially
        affected by edges effects.
    shift : dict
        A dictionary with two keys, 'x' and 'y'.  Key 'x' is an astropy
        quantities array of corresponding to the amount of shift in the
        x-direction (in arcseconds, assuming the helio-projective
        Cartesian co-ordinate system) that is applied to the input
        `~sunpy.map.MapCube`.  Key 'y' is an `~astropy.units.Quantity` array
        corresponding to the amount of shift in the y-direction (in arcseconds,
        assuming the helio-projective Cartesian co-ordinate system) that is
        applied to the input `~sunpy.map.MapCube`. The number of elements in
        each array must be the same as the number of maps in the
        `~sunpy.map.MapCube`.  If a shift is passed in to the function, that
        shift is applied to the input `~sunpy.map.MapCube` and the template
        matching algorithm is not used.

    The remaining keyword arguments are sent to `sunpy.image.coalignment.apply_shifts`.

    Returns
    -------
    output : `sunpy.map.MapCube`
        A `~sunpy.map.MapCube` that has co-aligned by matching the template.
    Examples
    --------
    >>> from sunpy.image.coalignment import mapcube_coalign_by_match_template as mc_coalign

    >>> coaligned_mc = mc_coalign(mc)   # doctest: +SKIP
    >>> coaligned_mc = mc_coalign(mc, layer_index=-1)   # doctest: +SKIP
    >>> coaligned_mc = mc_coalign(mc, clip=False)   # doctest: +SKIP
    >>> coaligned_mc = mc_coalign(mc, template=sunpy_map)   # doctest: +SKIP
    >>> coaligned_mc = mc_coalign(mc, template=two_dimensional_ndarray)   # doctest: +SKIP
    >>> coaligned_mc = mc_coalign(mc, func=np.log)   # doctest: +SKIP
    """
    return mapsequence_coalign_by_match_template(mc, template, layer_index,
                                                 func, clip, shift, **kwargs)


# Coalignment by matching a template
def mapsequence_coalign_by_match_template(mc, template=None, layer_index=0,
                                          func=_default_fmap_function, clip=True,
                                          shift=None, **kwargs):
    """
    Co-register the layers in a `~sunpy.map.MapSequence` according to a template
    taken from that `~sunpy.map.MapSequence`.  This method REQUIRES that
    scikit-image be installed. When using this functionality, it is a good idea
    to check that the shifts that were applied to were reasonable and expected.
    One way of checking this is to animate the original `~sunpy.map.MapSequence`,
    animate the coaligned `~sunpy.map.MapSequence`, and compare the differences you
    see to the calculated shifts.


    Parameters
    ----------
    mc : `sunpy.map.MapSequence`
        A `~sunpy.map.MapSequence` of shape (ny, nx, nt), where nt is the number of
        layers in the `~sunpy.map.MapSequence`.
    template : {None | sunpy.map.Map | `~numpy.ndarray`}
        The template used in the matching.  If an `~numpy.ndarray` is passed,
        the `~numpy.ndarray` has to have two dimensions.
    layer_index : int
        The template is assumed to refer to the map in the `~sunpy.map.MapSequence`
        indexed by the value of "layer_index".  Displacements of all maps in the
        `~sunpy.map.MapSequence` are assumed to be relative to this layer.  The
        displacements of the template relative to this layer are therefore
        (0, 0).
    func : function
        A function which is applied to the data values before the coalignment
        method is applied.  This can be useful in coalignment, because it is
        sometimes better to co-align on a function of the data rather than the
        data itself.  The calculated shifts are applied to the original data.
        Examples of useful functions to consider for EUV images are the
        logarithm or the square root.  The function is of the form
        func = F(data).  The default function ensures that the data are
        floats.
    clip : bool
        If True, then clip off x, y edges of the maps in the sequence that are
        potentially affected by edges effects.
    shift : dict
        A dictionary with two keys, 'x' and 'y'.  Key 'x' is an astropy
        quantities array of corresponding to the amount of shift in the
        x-direction (in arcseconds, assuming the helio-projective
        Cartesian co-ordinate system) that is applied to the input
        `~sunpy.map.MapSequence`.  Key 'y' is an `~astropy.units.Quantity` array
        corresponding to the amount of shift in the y-direction (in arcseconds,
        assuming the helio-projective Cartesian co-ordinate system) that is
        applied to the input `~sunpy.map.MapSequence`. The number of elements in
        each array must be the same as the number of maps in the
        `~sunpy.map.MapSequence`.  If a shift is passed in to the function, that
        shift is applied to the input `~sunpy.map.MapSequence` and the template
        matching algorithm is not used.

    The remaining keyword arguments are sent to `sunpy.image.coalignment.apply_shifts`.

    Returns
    -------
    output : `sunpy.map.MapSequence`
        A `~sunpy.map.MapSequence` that has co-aligned by matching the template.

    Examples
    --------
    >>> from sunpy.image.coalignment import mapsequence_coalign_by_match_template as mc_coalign

    >>> coaligned_mc = mc_coalign(mc)   # doctest: +SKIP
    >>> coaligned_mc = mc_coalign(mc, layer_index=-1)   # doctest: +SKIP
    >>> coaligned_mc = mc_coalign(mc, clip=False)   # doctest: +SKIP
    >>> coaligned_mc = mc_coalign(mc, template=sunpy_map)   # doctest: +SKIP
    >>> coaligned_mc = mc_coalign(mc, template=two_dimensional_ndarray)   # doctest: +SKIP
    >>> coaligned_mc = mc_coalign(mc, func=np.log)   # doctest: +SKIP
    """

    # Number of maps
    nt = len(mc.maps)

    # Storage for the pixel shifts and the shifts in arcseconds
    xshift_keep = np.zeros(nt) * u.pix
    yshift_keep = np.zeros_like(xshift_keep)

    if shift is None:
        shifts = calculate_match_template_shift(mc, template=template,
                                                layer_index=layer_index,
                                                func=func)
        xshift_arcseconds = shifts['x']
        yshift_arcseconds = shifts['y']
    else:
        xshift_arcseconds = shift['x']
        yshift_arcseconds = shift['y']

    # Calculate the pixel shifts
    for i, m in enumerate(mc):
        xshift_keep[i] = (xshift_arcseconds[i] / m.scale[0])
        yshift_keep[i] = (yshift_arcseconds[i] / m.scale[1])

    # Apply the shifts and return the coaligned mapsequence
    return apply_shifts(mc, -yshift_keep, -xshift_keep, clip=clip, **kwargs)
