import sunpy
import sunpy.map
import sunpy.cm
import sunpy.wcs

import numpy as np
import astropy.wcs
import astropy.io.fits as fits
import astropy.units as u

import scipy.ndimage as snd
import matplotlib.pyplot as plt
from glob import glob



def slit(in_files, xy1, xy2):
    """
    Returns an array with intensity along the slit on the y axis and time 
    along x.
    
    Parameters    
    ----------
    in_files : `sunpy.map.MapCube`
        A mapcube of the images you want to perform the slit analysis on. 
        Usually with x and y as space, and z as time.
    xy1 : List of x and y coordinates in `[x1, y1]`
        The x and y coordinates of the beginning of the slit. May be either 
        pixel coordinates or astropy units.
    xy2 : List of x and y coordinates
        The x and y coordinates of the end of the slit. May be either 
        pixel coordinates or astropy units.

        
    Returns
    -------
    out : numpy array
        A numpy array with intenisty in y and time in x.

    """
    
    # check the attributes of the coordinates
    if ((isinstance(xy1[0] and xy1[1], u.Quantity) and isinstance(xy2[0] and xy2[1], u.Quantity)) or
        (hasattr(xy1[0] and xy1[1], 'unit') and hasattr(xy2[0] and xy2[1], 'unit'))):
    
        if (xy1[0].unit.is_equivalent(in_files[0].units.x) and
            xy1[1].unit.is_equivalent(in_files[0].units.y)):
            units = 'data'
            
            # convert the world to pixel
            init_map = sunpy.map.Map(in_files[0])
            x1, y1 = init_map.data_to_pixel(xy1[0], xy1[1])
            x2, y2 = init_map.data_to_pixel(xy2[0], xy2[1])            
            
            
        elif xy1[0].unit.is_equivalent(u.pixel) and xy1[1].unit.is_equivalent(u.pixel):
            units = 'pixels'
        else:
            raise u.UnitsError("xy1 and xy2 must be "
                               "in units convertable to {} or {}".format(in_files[0].units['x'],
                                                                         u.pixel))
    else:
        raise TypeError("Arguments range_a and range_b to function submap "
                        "have an invalid unit attribute "
                        "You may want to pass in an astropy Quantity instead.")
        
    
    
    # call to the get pixel numbers routine
    slits = []
    slit = get_pixels_on_line(int(x1.value), int(y1.value), 
                              int(x2.value), int(y2.value))    
    


    for d_arr in in_files: 
        data = d_arr.data
        # get the initial slit pixels
        s_values = data[slit.T[0], slit.T[1]]
        # plus one
        s_p1_values = data[slit.T[0], slit.T[1]]
        # minus one
        s_m1_values = data[slit.T[0], slit.T[1]]
        # wap it all in one list
        slits.append([s_m1_values, s_values, s_p1_values])   
    
    end_array = np.array(slits)
    im_array = np.rot90(end_array)
    return(im_array)

def get_pixels_on_line(x1, y1, x2, y2, getvalues=True):
    """
    Uses Bresenham's line algorithm to enumerate the pixels along
    a line.
    (see http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm)
    If `getvalues`==False then it will return tuples of (x, y) coordinates
    instead of pixel values. This was taken from the package Ginga.
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




