import sunpy
import sunpy.map
import sunpy.cm
import sunpy.wcs

import numpy as np
import astropy.io.fits as fits
import astropy.units as u

import scipy.ndimage as snd
import matplotlib.pyplot as plt
from glob import glob

# import files command
files = glob('/storage2/SDO/jet/crop/304*.fits')
files.sort()


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
    
    slits = []
    slit = get_pixels_on_line(int(xy1[0].value), int(xy1[1].value), 
                              int(xy2[0].value), int(xy2[1].value))    
    


    for afile in files: 
        amap = sunpy.map.Map(afile)
        data = amap.data
        # get the initial slit pixels
        s_values = data[slit.T[0], slit.T[1]]
        # plus one
        s_p1_values = data[slit.T[0], slit.T[1]]
        # minus one
        s_m1_values = data[slit.T[0], slit.T[1]]
        # wap it all in one list
        slits.append([s_m1_values, s_values, s_p1_values])   
    
    end_array = np.array(slits)
    mean_array = end_array.mean(axis=1)
    array_maxes = mean_array.max(axis=0)
    norm_array = mean_array/array_maxes
    im_array = np.rot90(norm_array)
    return(im_array)

def get_pixels_on_line(x1, y1, x2, y2, getvalues=True):
    """Uses Bresenham's line algorithm to enumerate the pixels along
    a line.
    (see http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm)
    If `getvalues`==False then it will return tuples of (x, y) coordinates
    instead of pixel values.
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




