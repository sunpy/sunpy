"""
This module provides routines for the coalignment of images and mapcubes.

Currently this module provides image coalignment by template matching.
Which is partially inspired by the SSWIDL routine
`tr_get_disp.pro <http://hesperia.gsfc.nasa.gov/ssw/trace/idl/util/routines/tr_get_disp.pro>`_.

In this implementation, the template matching is handled via
the scikit-image routine :func:`skimage.feature.match_template`.

References
----------
Template matching algorithm:

 * http://scribblethink.org/Work/nvisionInterface/nip.html
 * J.P. Lewis, Fast Template Matching, Vision Interface 95, Canadian Image
   Processing and Pattern Recognition Society, Quebec City, Canada, May 15-19,
   1995, p. 120-123 http://www.scribblethink.org/Work/nvisionInterface/vi95_lewis.pdf.
"""

import numpy as np
from scipy.ndimage.interpolation import shift
from copy import deepcopy
from astropy import units as u

# Image co-registration by matching templates
from skimage.feature import match_template

# SunPy imports
from sunpy.map.mapbase import GenericMap

__author__ = 'J. Ireland'

__all__ = ['calculate_shift', 'clip_edges', 'calculate_clipping',
           'match_template_to_layer', 'find_best_match_location',
           'get_correlation_shifts', 'parabolic_turning_point',
           'repair_image_nonfinite', 'mapcube_coalign_by_match_template']


def _default_fmap_function(data):
    """This function ensures that the data are floats.  It is the default data
    manipulation function for the coalignment method.
    """
    return np.float64(data)


def calculate_shift(this_layer, template):
    """Calculates the pixel shift required to put the template in the "best"
    position on a layer.

    Parameters
    ----------
    this_layer : ndarray
        A numpy array of size (ny, nx), where the first two dimensions are
        spatial dimensions.

    template : ndarray
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
def clip_edges(data, yclips, xclips):
    """Clips off the y and x edges of a 2d array according to a list of pixel
    values.  This function is useful for removing data at the edge of
    2d images that may be affected by shifts from solar de-rotation and
    layer co-registration, leaving an image unaffected by edge effects.

    Parameters
    ----------
    data : ndarray
        A numpy array of shape (ny, nx).

    yclips : ndarray of astropy.units.Quantity
        The amount to clip in the y-direction of the data.

    xclips : ndarray of astropy.units.Quantity
        The amount to clip in the x-direction of the data.

    Returns
    -------
    image : ndarray
        A 2d image with edges clipped off according to the positive and
        negative ceiling values in the yclips and xclips arrays.
    """
    if not isinstance(yclips and xclips, u.Quantity):
        raise ValueError("Must be astropy QUantity")
    # Datacube shape
    ny = data.shape[0]
    nx = data.shape[1]
    return data[yclips[0].value: ny - yclips[1].value, xclips[0].value: nx - xclips[1].value]


#
# Return the upper and lower clipping values for the y and x directions an
# input set of pixel shifts y and x
#
def calculate_clipping(y, x):
    """Return the upper and lower clipping values for the y and x directions.

    Parameters
    ----------
    y : ndarray of type astropy.units.Quantity
        An array of pixel shifts in the y-direction for an image.

    x : ndarray of type astropy.units.Quantity
        An array of pixel shifts in the x-direction for an image.

    Returns
    -------
    clipping : ([int, int], [int, int]) of type astropy.Quantity
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
        columns).
    """
    if not isinstance(y, u.Quantity):
        raise ValueError("Must be astropy Quantites")
    if not isinstance(x, u.Quantity):
        raise ValueError("Must be astropy Quantites")
    return [_lower_clip(y.value), _upper_clip(y.value)] * u.pix, 
    [_lower_clip(x.value), _upper_clip(x.value)] * u.pix


#
# Helper functions for clipping edges
#
def _upper_clip(z):
    """Find smallest integer bigger than all the positive entries in the input
    array.
    """
    zupper = 0
    zcond = z >= 0
    if np.any(zcond):
        zupper = np.max(np.ceil(z[zcond]))
    return zupper


def _lower_clip(z):
    """Find smallest positive integer bigger than the absolute values of the
    negative entries in the input array.
    """
    zlower = 0
    zcond = z <= 0
    if np.any(zcond):
        zlower = np.max(np.ceil(-z[zcond]))
    return zlower


def match_template_to_layer(layer, template):
    """Calculate the correlation array that describes how well the template
    matches the layer. All inputs are assumed to be numpy arrays.  This
    function requires the "match_template" function in scikit image.

    Parameters
    ----------
    layer : ndarray
        A numpy array of size (ny, nx).

    template : ndarray
        A numpy array of size (N, M) where N < ny and M < nx.

    Returns
    -------
    correlationarray : ndarray
        A correlation array between the layer and the template.
        The values in the array range between 0 and 1.
    """
    return match_template(layer, template)


def find_best_match_location(corr):
    """Calculate an estimate of the location of the peak of the correlation
    result in image pixels.

    Parameters
    ----------
    corr : ndarray
        A 2-d correlation array.

    Returns
    -------
    shift : tuple of type astropy.Quantity
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
    """Estimate the location of the maximum of a fit to the input array.  The
    estimation in the x and y directions are done separately. The location
    estimates can be used to implement subpixel shifts between two different
    images.

    Parameters
    ----------
    array : ndarray
        An array with at least one dimension that has three elements.  The
        input array is at most a 3 x 3 array of correlation values calculated
        by matching a template to an image.

    Returns
    -------
    peakloc : tuple of type astropy.units.Quantity
        The (y, x) location of the peak of a parabolic fit, in image pixels.
    """
    # Check input shape
    ny = array.shape[0]
    nx = array.shape[1]
    if nx > 3 or ny > 3:
        print 'Input array is too big in at least one dimension. Returning Nones'
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
    """Find the location of the turning point for a parabola
    y(x) = ax^2 + bx + c, given input values y(-1), y(0), y(1).
    The maximum is located at x0 = -b / 2a .  Assumes
    that the input array represents an equally spaced sampling at the
    locations y(-1), y(0) and y(1).

    Parameters
    ----------
    y : ndarray
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
    """Return a new image in which all the nonfinite entries of the original
    image have been replaced by the local mean.

    Parameters
    ----------
    image : ndarray
        A two-dimensional ndarray.

    Returns
    -------
    repaired_image : ndarray
        A two-dimensional ndarray of the same shape as the input that has all
        the non-finite entries replaced by a local mean.  The algorithm
        repairs one non-finite entry at every pass.  At each pass, the next
        non-finite value is replaced by the mean of its finite valued nearest
        neighbours.
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


# Coalignment by matching a template
def mapcube_coalign_by_match_template(mc, template=None, layer_index=0,
                               func=_default_fmap_function, clip=True,
                               return_displacements_only=False,
                               apply_displacements=None,
                               with_displacements=False):
    """Co-register the layers in a mapcube according to a template taken from
    that mapcube.  This method REQUIRES that scikit-image be installed.
    When using this functionality, it is a good idea to check that the
    shifts that were applied to were reasonable and expected.  One way of
    checking this is to animate the original mapcube, animate the coaligned
    mapcube, and compare the differences you see to the calculated shifts.


    Parameters
    ----------
    mc : sunpy.map.MapCube
        A mapcube of shape (ny, nx, nt), where nt is the number of layers in
        the mapcube.

    template : {None | sunpy.map.Map | ndarray}
        The template used in the matching.  If an ndarray is passed, the
        ndarray has to have two dimensions.

    layer_index : int
        The template is assumed to refer to the map in the mapcube indexed by
        the value of "layer_index".  Displacements of all maps in the mapcube
        are assumed to be relative to this layer.  The displacements of the
        template relative to this layer are therefore (0, 0).

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
        If True, thenclip off x, y edges in the datacube that are potentially
        affected by edges effects.

    return_displacements_only : bool
        If True return ONLY the x and y displacements applied to the input
        data in units of arcseconds.  The return value is a dictionary of the
        form {"x": xdisplacement, "y": ydisplacement}.

    apply_displacements : {None | dict}
        If not None, then use the displacements supplied by the user.  Must be
        in the same format as that returned using the
        return_displacements_only option.  Can be used when you want to appl
        the same displacements to multiple mapcubes.

    with_displacements : bool
        If True, return the x and y displacements applied to the input data in
        the same format as that returned using the return_displacements_only
        option, along with the coaligned mapcube.  The format of the return is
        (mapcube, displacements).

    Returns
    -------
    output : {sunpy.map.MapCube | dict | tuple}
        The results of the mapcube coalignment.  The output depends on the
        value of the parameters "return_displacements_only" and
        "with_displacements".

    Examples
    --------
    >>> import numpy as np
    >>> from sunpy.image.coalignment import mapcube_coalign_by_match_template as mc_coalign
    >>> coaligned_mc = mc_coalign(mc)
    >>> coaligned_mc = mc_coalign(mc, layer_index=-1)
    >>> coaligned_mc = mc_coalign(mc, clip=False)
    >>> coaligned_mc = mc_coalign(mc, template=sunpy_map)
    >>> coaligned_mc = mc_coalign(mc, template=two_dimensional_ndarray)
    >>> coaligned_mc = mc_coalign(mc, func=np.log)
    >>> displacements = mc_coalign(mc, return_displacements_only=True)
    >>> coaligned_mc, displacements = mc_coalign(mc, with_displacements=True)
    >>> coaligned_mc = mc_coalign(mc, apply_displacements=displacements)
    """
    # Size of the data
    ny = mc.maps[layer_index].shape[0]
    nx = mc.maps[layer_index].shape[1]
    nt = len(mc.maps)

    # Storage for the pixel shifts and the shifts in arcseconds
    xshift_keep = np.zeros((nt))
    yshift_keep = np.zeros_like(xshift_keep)

    # Use the displacements supplied
    if apply_displacements is not None:
        xshift_arcseconds = apply_displacements["x"].value
        yshift_arcseconds = apply_displacements["y"].value
        for i, m in enumerate(mc.maps):
            xshift_keep[i] = xshift_arcseconds[i] / m.scale['x']
            yshift_keep[i] = yshift_arcseconds[i] / m.scale['y']
    else:
        xshift_arcseconds = np.zeros_like(xshift_keep)
        yshift_arcseconds = np.zeros_like(xshift_keep)

        # Calculate a template.  If no template is passed then define one
        # from the the index layer.
        if template is None:
            tplate = mc.maps[layer_index].data[ny / 4: 3 * ny / 4,
                                             nx / 4: 3 * nx / 4]
        elif isinstance(template, GenericMap):
            tplate = template.data
        elif isinstance(template, np.ndarray):
            tplate = template
        else:
            raise ValueError('Invalid template.')

        # Apply the function to the template
        tplate = func(tplate)

        # Match the template and calculate shifts
        for i, m in enumerate(mc.maps):
            # Get the next 2-d data array
            this_layer = func(m.data)

            # Calculate the y and x shifts in pixels
            yshift, xshift = calculate_shift(this_layer, tplate)

            # Keep shifts in pixels
            yshift_keep[i] = yshift.value
            xshift_keep[i] = xshift.value

        # Calculate shifts relative to the template layer
        yshift_keep = yshift_keep - yshift_keep[layer_index]
        xshift_keep = xshift_keep - xshift_keep[layer_index]

        for i, m in enumerate(mc.maps):
            # Calculate the shifts required in physical units, which are
            # presumed to be arcseconds.
            xshift_arcseconds[i] = xshift_keep[i] * m.scale['x']
            yshift_arcseconds[i] = yshift_keep[i] * m.scale['y']

    # Return only the displacements
    if return_displacements_only:
        return {"x": xshift_arcseconds * u.arcsec, "y": yshift_arcseconds * u.arcsec}

    # New mapcube for the new data
    newmc = deepcopy(mc)

    # Shift the data and construct the mapcube
    for i, m in enumerate(newmc.maps):
        shifted_data = shift(m.data, [-yshift_keep[i], -xshift_keep[i]])
        if clip:
            yclips, xclips = calculate_clipping(yshift_keep, xshift_keep)
            shifted_data = clip_edges(shifted_data, yclips.value, xclips.value)

        # Update the mapcube image data
        newmc.maps[i].data = shifted_data

        # Adjust the positioning information accordingly.
        newmc.maps[i].meta['crpix1'] = newmc.maps[i].meta['crpix1'] + xshift_arcseconds[i]
        newmc.maps[i].meta['crpix2'] = newmc.maps[i].meta['crpix2'] + yshift_arcseconds[i]

    # Return the mapcube, or optionally, the mapcube and the displacements
    # used to create the mapcube.
    if with_displacements:
        return newmc, {"x": xshift_arcseconds * u.arcsec, "y": yshift_arcseconds * u.arcsec}
    else:
        return newmc

