import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

def darklimb_u(wavelength):
    """
    Returns the constant, u, for a given wavelength for input into the limb 
    darkening function.

    Parameters
    ----------
    wavelength : `~sunpy.map.GenericMap.wavelength`
        The wavelength for which the limb darkening is to be calculated.

    Notes
    -----
    Coefficients for 5th order polynomial fit to limb darkening
    constant `u` was done on data from Astrophysical
    Quantities by Allen.

    Limitations
    -----------
    Only applies to wavelengths in the range 4000<lambda<15000 Angstrom

    Reference
    ---------
    Cox, Arthur N., Allen's astrophysical quantities, 4th ed. Publisher:
        New York: AIP Press; Springer, 2000. ISBN: 038798746

    """

    au = -8.9829751
    bu = 0.0069093916
    cu = -1.8144591e-6
    du = 2.2540875e-10
    eu = -1.3389747e-14
    fu = 3.0453572e-19

    a=np.array((au,bu,cu,du,eu,fu))

    wavelength = wavelength.to_value()
    wavelength_series = np.array([wavelength ** n for n in range(len(a))])

    ul = np.sum(a * wavelength_series)

    return ul

def darklimb_v(wavelength):
    """
    Returns the constant, v, for a given wavelength for input into the limb 
    darkening function.

    Parameters
    ----------
    wavelength : `~sunpy.map.GenericMap.wavelength`
        The wavelength for which the limb darkening is to be calculated.

    Notes
    -----
    Coefficients for 5th order polynomial fit to limb darkening
    constant `v` was done on data from Astrophysical
    Quantities by Allen.

    Limitations
    -----------
    Only applies to wavelengths in the range 4000<lambda<15000 Angstrom

    Reference
    ---------
    Cox, Arthur N., Allen's astrophysical quantities, 4th ed. Publisher:
        New York: AIP Press; Springer, 2000. ISBN: 038798746

    """

    av = 9.2891180
    bv = -0.0062212632
    cv = 1.5788029e-6
    dv = -1.9359644e-10
    ev = 1.1444469e-14
    fv = -2.599494e-19

    a=np.array((av,bv,cv,dv,ev,fv))

    wavelength = wavelength.to_value()
    wavelength_series = np.array([wavelength ** n for n in range(len(a))])

    vl = np.sum(a * wavelength_series)

    return vl

def get_limb_dark(sun_map):
    """
    Calculates limb darkening function for a given sun_map

    Parameters
    ----------
    sun_map : `~sunpy.map.GenericMap`
        The sun map for which the limb darkening is to be calculated.

    Limitations
    -----------
    Only applies to wavelengths in the range 4000<lambda<15000 Angstrom

    Returns
    -------
    limbfilt: `~numpy.array`
        Numpy array representation of limb darkening function.

    Notes
    -----
    IDL code equivalent:
    http://hesperia.gsfc.nasa.gov/ssw/packages/nrl/idl/nrlgen/analysis/limb_dark.pro

    Reference
    ---------
    Cox, Arthur N., Allen's astrophysical quantities, 4th ed. Publisher:
        New York: AIP Press; Springer, 2000. ISBN: 038798746

    """

    if not (4000 * u.AA < sun_map.wavelength < 15000 * u.AA):
        raise ValueError("The wavelength of the observation must be \
                          between 4000 and 15000 Angstrom.")

    ul = darklimb_u(sun_map.wavelength)
    vl = darklimb_v(sun_map.wavelength)

    x_center, y_center = sun_map.world_to_pixel(SkyCoord(0, 0, unit=u.arcsec, frame=sun_map.coordinate_frame))
    x_center = x_center.to(u.pixel)
    y_center = y_center.to(u.pixel)
    x_dim = sun_map.dimensions.x.to(u.pixel).value
    y_dim = sun_map.dimensions.y.to(u.pixel).value
    radius = (sun_map.rsun_obs.to(u.arcsec).value /
              sun_map.scale.axis1.to(u.arcsec / u.pixel).value)

    x_grid, y_grid = np.meshgrid(np.arange(int(x_dim)), np.arange(int(y_dim)),sparse=True)

    x_2 = (x_grid - x_center.value) ** 2
    y_2 = (y_grid - y_center.value) ** 2

    dist_grid = np.sqrt(x_2 + y_2) / radius
    dist_grid[np.where(dist_grid>1)] = 0

    limbfilt = 1 - ul - vl + ul* np.cos(np.arcsin(dist_grid)) + vl* np.cos(np.arcsin(dist_grid))**2

    return limbfilt

def correct_limb_darkening(sun_map):
    """
    Removes limb darkening effects from input sun_map.

    Parameters
    ----------
    sun_map : `~sunpy.GenericMap`
        The sun map from which the limb darkening effects are to be removed.

    Returns
    -------
    new_map: `~sunpy.GenericMap`
        New instance of input sun map with limb darkening effects removed.

    Limitations
    -----------
    Only applies to wavelengths in the range 4000<lambda<15000 Angstrom

    Examples
    --------
    >>> import sunpy.map
    >>> from sunpy.phyisics import limb_darkening
    >>> mymap = sunpy.map.Map('test.fits')
    >>> newmap = limb_darkening.correct_limb_darkening(mymap)
    >>> mymap.peek()
    >>> newmap.peek()
    """
    new_data = sun_map.data / get_limb_dark(sun_map)
    new_meta = sun_map.meta
    new_meta['history'] += "Limb Darkening effects removed. The limb darkening function uses a 5th order polynomial fitting\
		                    to the limb darkening constants obtained from Astrophysical Quantities."

    return (sun_map._new_instance(new_data, new_meta, sun_map.plot_settings))
