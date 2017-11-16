import copy as c
import sunpy.map
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

__author__ = ["Norbert Gyenge"]


def limbdark(observation, limb_cut=True):
    """
    Calculates the limb darkening function for a specified wavelength
    and removes any limb darkening effect from the input image.
    A 5th order polynomial is fit to the limb darkening constants obtained from
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
        observation - Corrected sunpy map object
        The original map will be kept.
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
        >>> import sunpy.map
        >>> mymap = sunpy.map.Map('test.fits')
        >>> newmap = limb_darkening(mymap, limb_cut=True)
        >>> mymap.peek()
        >>> newmap.peek()"""

    if observation.wavelength.value < 4000 or observation.wavelength.value > 15000:
        raise ValueError("The wavelength of the observation must be between 4000 and 15000 A.")

    observation = c.copy(observation)
    wavelength = observation.wavelength.value
    x_center, y_center = observation.world_to_pixel(
        SkyCoord(0, 0, unit=u.arcsec, frame=observation.coordinate_frame))
    x_dim, y_dim = observation.dimensions.x.value, observation.dimensions.y.value
    radius = observation.rsun_obs.value / observation.scale.axis1.value

    a = np.array([-8.9829751, 0.0069093916, -1.8144591e-6, 2.2540875e-10,
                  -1.3389747e-14, 3.0453572e-19])
    b = np.array([9.2891180, -0.0062212632, 1.5788029e-6, -1.9359644e-10,
                  1.1444469e-14, -2.599494e-19])

    wavelength = [1, wavelength, wavelength ** 2, wavelength ** 3,
                  wavelength ** 4, wavelength ** 5]

    ul = sum(a * wavelength)
    vl = sum(b * wavelength)

    x_grid, y_grid = np.meshgrid(np.arange(int(x_dim)), np.arange(int(y_dim)))
    x_2 = (x_grid - x_center.value) ** 2
    y_2 = (y_grid - y_center.value) ** 2

    dist_grid = np.sqrt(x_2 + y_2)
    dist_grid = dist_grid / radius

    e1 = 1 - ul - vl + ul * np.cos(np.arcsin(dist_grid))
    e2 = vl * (np.cos(np.arcsin(dist_grid)) ** 2)
    limbfilt = e1 + e2

    res = observation.data / limbfilt
    new_observation = sunpy.map.Map(res.astype(np.float32), observation.meta)
    new_observation.meta['history'] += 'DARKLIMB CORRECTED'

    if limb_cut is True:
        new_observation.data[dist_grid > radius] = np.nan

    return new_observation
