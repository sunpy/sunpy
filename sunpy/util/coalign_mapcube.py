import numpy as np
from sunpy.map import Map

# For faster image co-registration
from skimage.feature import match_template

# Shift an image by a given amount - subpixel shifts are permitted
from scipy.ndimage.interpolation import shift


def default_data_manipulation_function(data):
    """
    This function ensures that the data are floats.  It is the default data
    manipulation function for the coalignment method.
    """
    return 1.0 * data


#
# Coalign a mapcube
#
def coalign_mapcube(mc,
                    layer_index=0,
                    func=default_data_manipulation_function,
                    clip=False):
    """
    Co-register the layers in a mapcube according to a template taken from
    that mapcube.

    Input
    -----
    mc : a mapcube of shape (ny, nx, nt), where nt is the number of
         layers in the mapcube.

    layer_index : the layer in the mapcube from which the template will be
                  extracted.

    func: a function which is applied to the data values before the
          coalignment method is applied.  This can be useful in coalignment,
          because it is sometimes better to co-align on a function of the data
          rather than the data itself.  The calculated shifts are applied to
          the original data.  Useful functions to consider are the log of the
          image data, or 1 / data. The function is of the form func = F(data).
          The default function ensures that the data are floats.

    clip : clip off x, y edges in the datacube that are potentially affected
            by edges effects.

    Output
    ------
    datacube : the input datacube each layer having been co-registered against
               the template.

    y_shift_keep : a one dimensional array of length nt with the pixel
                     y-displacements relative position of the template at the
                     value layer_index.  Note that y_displacement[layer_index]
                     is zero by definition.

    x_shift_keep: a one dimensional array of length nt with the pixel
                     x-displacements relative position of the template at the
                     value layer_index.  Note that x_displacement[layer_index]
                     is zero by definition.
    """
    # Size of the data
    ny = mc.maps[layer_index].shape[0]
    nx = mc.maps[layer_index].shape[1]
    nt = len(mc.maps)

    # Storage for the shifted data and the pixel shifts
    shifted_datacube = np.zeros((ny, nx, nt))
    xshift_keep = np.zeros((nt))
    yshift_keep = np.zeros((nt))

    # Calculate a template
    template = func(mc.maps[layer_index].data[ny / 4: 3 * ny / 4,
                                         nx / 4: 3 * nx / 4])

    for i, m in enumerate(mc.maps):
        # Get the next 2-d data array
        this_layer = func(m.data)

        # Calculate the y and x shifts in pixels
        yshift, xshift = calculate_shift(this_layer, template)

        # Keep shifts in pixels
        yshift_keep[i] = yshift
        xshift_keep[i] = xshift

    # Calculate shifts relative to the template layer
    yshift_keep = yshift_keep - yshift_keep[layer_index]
    xshift_keep = xshift_keep - xshift_keep[layer_index]

    # Shift the data
    for i, m in enumerate(mc.maps):
        shifted_datacube[:, :, i] = shift(m.data, [-yshift_keep[i], -xshift_keep[i]])

    # Clip the data if requested
    if clip:
        shifted_datacube = clip_edges(shifted_datacube, yshift_keep, xshift_keep)

    # Create a new mapcube.  Adjust the positioning information accordingly.
    new_cube = []
    for i, m in enumerate(mc.maps):
        print i
        new_meta = m.meta.copy()
        new_meta['xcen'] = new_meta['xcen'] + xshift_keep[i] * m.scale['x']
        new_meta['ycen'] = new_meta['ycen'] + yshift_keep[i] * m.scale['y']
        print new_meta['xcen'], new_meta['ycen']
        new_map = Map(shifted_datacube[:, :, i], new_meta)

        # Store the new map in a list
        new_cube.append(new_map)

    # Return the cube and the pixel displacements
    return Map(new_cube, cube=True), yshift_keep, xshift_keep


def calculate_shift(this_layer, template):
    """
    Calculates the pixel shift required to put the template in the "best"
    position on a layer.

    Inputs
    ------
    template : a numpy array of size (N, M) where N < ny and M < nx .

    this_layer : a numpy array of size (ny, nx), where the first two
               dimensions are spatial dimensions.

    Outputs
    -------
    yshift, xshift : pixel shifts relative to the offset of the template to
                     the input array.
    """
    # Repair any NANs, Infs, etc in the layer and the template
    this_layer = repair_nonfinite(this_layer)
    template = repair_nonfinite(template)

    # Calculate the correlation array matching the template to this layer
    corr = match_template_to_layer(this_layer, template)

    # Calculate the y and x shifts in pixels
    return find_best_match_location(corr)


#
# Remove the edges of a datacube
#
def clip_edges(datacube, y, x):
    """
    Clips off the y and x edges of a datacube according to a list of pixel
    values.  Positive pixel values will clip off the datacube at the upper end
    of the range.  Negative values will clip off values at the lower end of
    the  range.  This function is useful for removing data at the edge of
    datacubes that may be affected by shifts from solar de-rotation and
    layer co-registration, leaving a datacube unaffected by edge effects.

    Input
    -----
    datacube : a numpy array of shape (ny, nx, nt), where nt is the number of
               layers in the datacube.

    y : a numpy array of pixel values that correspond to how much to pixel
        clip values in the x-direction.

    x : a numpy array of pixel values that correspond to how much to pixel
        clip values in the y-direction.

    Output
    ------
    A datacube with edges clipd off according to the positive and negative
    ceiling values in the y and x arrays.
    """

    # Datacube shape
    ny = datacube.shape[0]
    nx = datacube.shape[1]
    nt = datacube.shape[2]

    # maximum to clip off in the x-direction
    xupper = _upper_clip(x)
    xlower = _lower_clip(x)

    # maximum to clip off in the y direction
    yupper = _upper_clip(y)
    ylower = _lower_clip(y)

    return datacube[ylower: ny - yupper - 1, xlower: nx - xupper - 1, 0: nt]


#
# Helper functions for clipping edges
#
def _upper_clip(z):
    zupper = 0
    zcond = z >= 0
    if np.any(zcond):
        zupper = np.max(np.ceil(z[zcond]))
    return zupper


def _lower_clip(z):
    zlower = 0
    zcond = z <= 0
    if np.any(zcond):
        zlower = np.max(np.ceil(-z[zcond]))
    return zlower


def match_template_to_layer(layer, template):
    """
    Calculate the correlation array that describes how well the template
    matches the layer.
    All inputs are assumed to be numpy arrays.

    Inputs
    ------
    template : a numpy array of size (N, M) where N < ny and M < nx .

    layer : a numpy array of size (ny, nx), where the first two
               dimensions are spatial dimensions.

    Outputs
    -------
    A cross-correlation array.  The values in the array range between 0 and 1.

    Requires
    --------
    This function requires the "match_template" function in scikit image.

    """
    return match_template(layer, template)


def find_best_match_location(corr):
    """
    Calculate an estimate of the location of the peak of the correlation
    result.

    Inputs
    ------
    corr : a 2-d correlation array.

    Output
    ------
    y, x : the shift amounts.  Subpixel values are possible.

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
    y_shift_relative_to_correlation_array = y_shift_relative_to_maximum + cor_max_y
    x_shift_relative_to_correlation_array = x_shift_relative_to_maximum + cor_max_x

    return y_shift_relative_to_correlation_array, x_shift_relative_to_correlation_array


def get_correlation_shifts(array):
    """
    Estimate the location of the maximum of a fit to the input array.  The
    estimation in the x and y directions are done separately. The location
    estimates can be used to implement subpixel shifts between two different
    images.

    Inputs
    ------
    array : an array with at least one dimension that has three elements.  The
            input array is at most a 3 x 3 array of correlation values
            calculated by matching a template to an image.

    Outputs
    -------
    y, x : the location of the peak of a parabolic fit.

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

    return y_location, x_location


def parabolic_turning_point(y):
    """
    Find the location of the turning point for a parabola f(x) = ax^2 + bx + c
    The maximum is located at x0 = -b / 2a .  Assumes that the input array
    represents an equally spaced sampling at the locations f(-1), f(0) and
    f(1).

    Input
    -----
    An one dimensional numpy array of shape 3

    Output
    ------
    A digit, the location of the parabola maximum.

    """
    numerator = -0.5 * y.dot([-1, 0, 1])
    denominator = y.dot([1, -2, 1])
    return numerator / denominator


def repair_nonfinite(z):
    """
    Replace all the nonfinite entries in a layer with the local mean.  There is
    probably a much smarter way of doing this.
    """
    nx = z.shape[1]
    ny = z.shape[0]
    bad_index = np.where(np.logical_not(np.isfinite(z)))
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
        subarray = z[y - 1: y + 2, x - 1: x + 2]
        z[by, bx] = np.nanmean(subarray * np.isfinite(subarray))
        bad_index = np.where(np.logical_not(np.isfinite(z)))
    return z


#
# Test functions for the coalignment code
#
def test_parabolic_turning_point():
    pass


def test_repair_nonfinite():
    pass


def test_get_correlation_shifts():
    pass