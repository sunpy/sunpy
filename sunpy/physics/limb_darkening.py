import sunpy.map
import numpy as np
from astropy import units as u
import copy as c
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord

# mymap = sunpy.map.Map('./*.fits')
# m = mymap[0].rotate(angle = mymap[0].meta.get('crota2') * u.deg)

def limbdark(m):
    m = c.copy(m)
    wavelength = m.wavelength.value
    x_center, y_center = m.world_to_pixel(SkyCoord(0, 0, unit=u.arcsec, frame=m.coordinate_frame))
    x_dim, y_dim = m.dimensions.x.value, m.dimensions.y.value
    radius = m.rsun_obs.value / m.scale.axis1.value
    
    a = np.array([-8.9829751, 0.0069093916, -1.8144591e-6, 2.2540875e-10,
                  -1.3389747e-14, 3.0453572e-19])
    b = np.array([9.2891180, -0.0062212632, 1.5788029e-6, -1.9359644e-10,
                  1.1444469e-14, -2.599494e-19])
    
    wavelength = [1, wavelength, wavelength ** 2, wavelength ** 3,
                  wavelength ** 4, wavelength ** 5]
    
    ul = sum(a * wavelength)
    vl = sum(b * wavelength)
    
    #x_grid, y_grid = np.mgrid[0:int(x_dim), 0:int(y_dim)]
    x_grid, y_grid = np.meshgrid(np.arange(int(x_dim)), np.arange(int(y_dim)))
    x_2 = (x_grid - x_center.value) ** 2
    y_2 = (y_grid - y_center.value) ** 2
    
    dist_grid = np.sqrt(x_2 + y_2)
    dist_grid = dist_grid / radius
    
    e1 = 1 - ul - vl + ul * np.cos(np.arcsin(dist_grid))
    e2 = vl * (np.cos(np.arcsin(dist_grid)) ** 2)
    limbfilt = e1 + e2
    
    rest = m.data / limbfilt
    mld = sunpy.map.Map(rest[1:-1,1:-1].astype(np.float32), m.meta)
    mld.meta['history'] += 'DARKLIMB CORRECTED'
    
    return mld
# plt.ion()
# plt.figure(2)
# plt.subplot(121)
# plt.subplot(121, projection=m)
# m.plot()
# ax.set_autoscale_on(False)
# plt.subplot(122, projection=m)
# ax.set_autoscale_on(False)
# plt.imshow(rest, cmap='gray')

#m.data = m.data / limbfilt

#if limb_cut is True:
#    m.data[dist_grid > radius] = np.nan
