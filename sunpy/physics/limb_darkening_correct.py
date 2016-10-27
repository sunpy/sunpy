import numpy as np
from astropy import units as u
__author__ = ["Norbert Gyenge"]


def limb_darkening_correct(observation, limb_cut=True):
    """Calculate limb darkening function for specified wavelength
    and use to remove limb darkening effects from input image.
    The limb darkening function uses a 5th ordern polynomial
    fitting to the limb darkening constants obtained from
    Astrophysical Quantities by Allan. Only in the wavelength
    range between 4000 and 15000 A.

    Parameters
    ----------
        observation - SunPy Map object
            SDO observation

        limb_cut - Default: True
            Background (r>1) values will be repleaced by np.nan

    Returns
    -------
        observation - Corrected map

    Notes
    -----
        IDL code equavalent:
        http://hesperia.gsfc.nasa.gov/ssw/packages/nrl/idl/nrlgen/analysis/limb_dark.pro

    Reference
    ---------
        Cox, Arthur N., Allen's astrophysical quantities, 4th ed. Publisher:
            New York: AIP Press; Springer, 2000. ISBN: 038798746

    Examples
    --------
        >>> mymap = sunpy.map.Map('test.fits')
        >>> limb_darkening_correct(mymap, limb_cut=True)"""

    wavelength = observation.wavelength.value

    if wavelength < 4000 or wavelength > 15000:
        raise ValueError("The wavelength of the observation must be between 4000 and 15000 A.")

    x_center, y_center = observation.data_to_pixel(0*u.arcsec, 0*u.arcsec)
    radius = observation.rsun_obs / observation.scale[0]
    x_dim, y_dim = observation.dimensions

    x_2 = np.power(np.arange(x_dim.value) - x_center.value, 2)
    y_2 = np.power(np.arange(y_dim.value) - y_center.value, 2)

    a = np.array([-8.9829751, 0.0069093916, -1.8144591e-6, 2.2540875e-10,
                  -1.3389747e-14, 3.0453572e-19])
    b = np.array([9.2891180, -0.0062212632, 1.5788029e-6, -1.9359644e-10,
                  1.1444469e-14, -2.599494e-19])

    wavelength = [1, wavelength, np.power(wavelength, 2),np.power(wavelength, 3),
                  np.power(wavelength, 4), np.power(wavelength, 5)]

    ul = sum(a * wavelength)
    vl = sum(b * wavelength)

    dist_grid = np.empty((int(x_dim.value), int(y_dim.value)), dtype=np.float)
    
    for i in range(0, int(y_dim.value)-1):
        dist_grid[:][i] = np.sqrt(y_2[i] + x_2)
        
    dist_grid = dist_grid / radius.value

    e1 = 1 - ul - vl + ul * np.cos(np.arcsin(dist_grid))
    e2 = vl * np.power(np.cos(np.arcsin(dist_grid)), 2)
    limbfilt = e1 + e2

    observation.data = observation.data / limbfilt

    if limb_cut is True:
        observation.data[dist_grid >= 1] = np.nan

    return observation
