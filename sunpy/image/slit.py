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



# define the range to preform the slit on in pixels
amap = sunpy.map.Map(files[100])
bl = amap.data_to_pixel(908.0*u.arcsec, -245.0*u.arcsec)
tr = amap.data_to_pixel(934.5*u.arcsec, -235.5*u.arcsec)

x_pixel = range(int(bl[0].value), int(tr[0].value))
ind_len = len(x_pixel)
y_pixel = np.linspace(int(bl[1].value), int(tr[1].value), num=ind_len)

inds = np.array([x_pixel,y_pixel])


# calculate the intensity along each slit
# append these to a list
slits = []
slit = get_pixels_on_line(int(bl[0].value), int(bl[1].value), int(tr[0].value), int(tr[1].value))
slit_p1 = slit + [-1,+1]
slit_m1 = slit + [+1,-1]
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





the_map = sunpy.map.Map(files[140])
plt.figure()
the_map.plot()
plt.plot(slit[:,0], slit[:,1])
plt.plot(slit_p1[:,0], slit_p1[:,1])
plt.plot(slit_m1[:,0], slit_m1[:,1])
plt.show()



plt.figure()
plt.imshow(im_array, origin='lower', cmap='sdoaia304')
plt.xlabel('Time')
plt.ylabel('Slit Intensity')
#plt.savefig('/data/SDO/images/211_slit_2.png')
plt.show()

