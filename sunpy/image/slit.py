import numpy as np

import astropy.units as u


def slit(mcube_in, range_a, range_b, N_slits=0, shift_x=+1, shift_y=-1,):
    """
    Returns an array with intensity along the slit on the y axis and time
    along x.

    Parameters
    ----------
    mcube_in : `sunpy.map.MapCube`
        A mapcube of the images you want to perform the slit analysis on.
        Usually with x and y as space, and z as time.
        The x and y axes must be equal for all instances within the mapcube.
    range_a : `astropy.units.Quantity`
        A list of two `astropy.unit.Quantity` objects representing x1 and x2,
        start and end of slit in x.
    range_b : `astropy.units.Quantity`
        A list of two `astropy.unit.Quantity` objects representing x2 and y2,
        start and end of slit in y.
    shift_x : `int`
        Possible values, `-1`, `0`, `+1`
        The displacement from the origin on the slit in x. These extra slits
        calcuate a mean along each row to make the positioning of the original
        slit less sensetive. Both `x` and `y` cannot be `0`. Default = `+1`
    shift_y : `int`
        Possible values, `-1`, `0`, `+1`
        The displacement from the origin in `y`. Default = `-1`
    N_slits : `int`
        Number of deviations from central slit, e.g. to use 5 slits `N = 2`.
        Default = 0

    Returns
    -------
    im_array : `numpy.ndarray`
        A numpy array with intenisty in y and time in x.
    slit : `numpy.ndarray`
        A numpy array representing the indicies of the central slit.

    """

    # check the attributes of the coordinates
    if ((isinstance(range_a and range_b, u.Quantity) or
        (hasattr(range_a and range_b, 'unit')))):
        if (range_a.unit.is_equivalent(mcube_in[0].units.x) and
            range_b.unit.is_equivalent(mcube_in[0].units.y)):
            # convert the world to pixel
            init_map = mcube_in[0]
            c_x1, c_y1 = init_map.data_to_pixel(range_a[0], range_b[0])
            c_x2, c_y2 = init_map.data_to_pixel(range_a[1], range_b[1])

        elif range_a.unit.is_equivalent(u.pixel) and range_b.unit.is_equivalent(u.pixel):
            c_x1, c_y1 = range_a[0], range_b[0]
            c_x2, c_y2 = range_a[1], range_b[1]
        else:
            raise u.UnitsError("xy1 and xy2 must be "
                               "in units convertable to {} or {}".format(mcube_in[0].units['x'],
                                                                         u.pixel))
    else:
        raise TypeError("Arguments range_a and range_b to function submap "
                        "have an invalid unit attribute "
                        "You may want to pass in an astropy Quantity instead.")

    if (shift_x and shift_y) == 0:
        raise ValueError("`shift_x` and `shift_y` are both equal to zero, therefore"
                         " not shift is possible")


    # extract the intensities from the data using function
    intensity_inds = slit_count(c_x1, c_y1, c_x2, c_y2, shift_y, shift_x, N_slits)

    # create an array of the data to index with the resualts of the slit file
    all_data = mcube_in.as_array()

    del_s_x = np.array([intensity_inds[:,:,0]])
    del_s_y = np.array([intensity_inds[:,:,1]])

    im_arr = all_data[del_s_x, del_s_y]
    im_out = im_arr.mean(axis=1)
    return im_out, intensity_inds[0]


def slit_count(bx, by, tx, ty, shift_x, shift_y, N):
    """
    Returns a list of lists associated with increments away from the inital
    slit. It transorms the slit array

    Parameters
    ----------
    bx, by, tx, ty : `astropy.unit.pixel`
        X and Y coordinates of the base and tip of the slit.
    shift_x : `int`
        Possible values, `-1`, `0`, `+1`
    shift_y : `int`
        Possible values, `-1`, `0`, `+1`
    N : `int`
        Number of slits away from central slit.

    Returns
    -------
    ind_array : `numpy.ndarray`
        Array of indicies representing the indicies of the interated slits
    """

    # call to the get pixel numbers routine
    if isinstance([bx, by, tx, ty], u.Unit) == True:
        init_slit = get_pixels_on_line(int(bx.value), int(by.value),
                                       int(tx.value), int(ty.value))
    else:
        init_slit = get_pixels_on_line(bx, by, tx, ty)

    shift_ind = np.array([shift_y, shift_x])
    shift_inds = [np.array([0, 0])]
    for i in range(1, N+1):
        shift_inds.append(shift_ind*i)
        shift_inds.append(shift_ind*(-i))

    list_slits = []
    for i in range(len(shift_inds)):
        slit_int = init_slit + shift_inds[i]
        list_slits.append(slit_int)

    return np.array(list_slits)


def get_pixels_on_line(x1, y1, x2, y2):
    """
    Returns an array of all pixel coordinates which the line defined by `x1, y1` and
    `x2, y2` crosses. Uses Bresenham's line algorithm to enumerate the pixels along
    a line. This was adapted from ginga

    Parameters
    ----------
    x1, y1, x2, y2 :`int`

    References
    ----------
    | https://github.com/ejeschke/ginga/blob/master/ginga/BaseImage.py#L387
    | http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
    | https://ginga.readthedocs.org/en/latest/

    """
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    if x1 < x2:
        sx = 1
    else:
        sx = -1
    if y1 < y2:
        sy = 1
    else:
        sy = -1
    err = dx - dy

    res = []
    x, y = x1, y1
    while True:
        res.append((x, y))
        if (x == x2) and (y == y2):
            break
        e2 = 2 * err
        if e2 > -dy:
            err = err - dy
            x += sx
        if e2 <  dx:
            err = err + dx
            y += sy

    return np.array(res)
