import sunpy.map

import numpy as np

import astropy.units as u


def slit(mcube_in, range_a, range_b):
    """
    Returns an array with intensity along the slit on the y axis and time
    along x.

    Parameters
    ----------
    mcube_in : `sunpy.map.MapCube`
        A mapcube of the images you want to perform the slit analysis on.
        Usually with x and y as space, and z as time.
        Within the mapcube x and y should be uniformly equal
    range_a : `astropy.units.Quantity`
        A list of two `astropy.unit.Quantity` objects representing x1 and x2,
        start and end of slit in x.
    range_b : `astropy.units.Quantity`
        A list of two `astropy.unit.Quantity` objects representing x2 and y2,
        start and end of slit in y.

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



    # call to the get pixel numbers routine
    slit = get_pixels_on_line(int(c_x1.value), int(c_y1.value),
                              int(c_x2.value), int(c_y2.value))
    # plus one, minus one, to get an average
    slit_p1 = slit + [-1,+1]
    slit_m1 = slit + [+1,-1]
    

    
    
    
    
    
    
    for d_arr in mcube_in:
        data = d_arr.data
        # get the initial slit pixels
        s_values = data[slit.T[0], slit.T[1]]
        # plus one
        s_p1_values = data[slit_p1.T[0], slit_p1.T[1]]
        # minus one
        s_m1_values = data[slit_m1.T[0], slit_m1.T[1]]
        # wap it all in one list
        slit_I.append([s_m1_values, s_values, s_p1_values])

    end_array = np.array(slit_I)
    mean_arr = end_array.mean(axis=1)
    im_array = np.rot90(mean_arr)
    return([im_array, slit])



def get_pixels_on_line(x1, y1, x2, y2):
    """
    Returns an array of all pixel coordinates which the line defined by `x1, y1` and
    `x2, y2` crosses. Uses Bresenham's line algorithm to enumerate the pixels along
    a line. This was adapted from ginga

    Parameters
    ----------
    x1, y1, x2, y2 : `<type 'int'>`

    If `getvalues`==False then it will return tuples of (x, y) coordinates
    instead of pixel values.

    References
    ----------
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