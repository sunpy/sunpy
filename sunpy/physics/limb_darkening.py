import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import copy


@u.quantity_input
def _darklimb_u(wavelength: u.AA):
    """
    Returns the constant, u, for a given wavelength for input into the limb
    darkening function.

    Coefficients for 5th order polynomial fit to limb darkening
    constant ``u`` was done on data from Astrophysical Quantities by Allen.

    Only applies to wavelengths in the range 4000<lambda<15000 Angstrom.

    Parameters
    ----------
    wavelength : `~astropy.units.Quantity`
        The wavelength for which the limb darkening is to be calculated.

    References
    ----------
    Cox, Arthur N., Allen's astrophysical quantities, 4th ed. Publisher:
        New York: AIP Press; Springer, 2000. ISBN: 038798746
    DOI : https://doi.org/10.1063/1.1325201
    """

    au = -8.9829751
    bu = 0.0069093916
    cu = -1.8144591e-6
    du = 2.2540875e-10
    eu = -1.3389747e-14
    fu = 3.0453572e-19

    a = np.array((au, bu, cu, du, eu, fu))

    wavelength_series = np.array([wavelength.value ** n for n in range(len(a))])

    ul = np.sum(a * wavelength_series)

    return ul

@u.quantity_input
def _darklimb_v(wavelength: u.AA):
    """
    Returns the constant, v, for a given wavelength for input into the limb
    darkening function.

    Coefficients for 5th order polynomial fit to limb darkening
    constant `v` was done on data from Astrophysical Quantities by Allen.

    Only applies to wavelengths in the range 4000<lambda<15000 Angstrom.

    Parameters
    ----------
    wavelength : `~astropy.units.Quantity`
        The wavelength for which the limb darkening is to be calculated.

    References
    ----------
    Cox, Arthur N., Allen's astrophysical quantities, 4th ed. Publisher:
        New York: AIP Press; Springer, 2000. ISBN: 038798746
    DOI : https://doi.org/10.1063/1.1325201
    """

    av = 9.2891180
    bv = -0.0062212632
    cv = 1.5788029e-6
    dv = -1.9359644e-10
    ev = 1.1444469e-14
    fv = -2.599494e-19

    a = np.array((av, bv, cv, dv, ev, fv))

    wavelength_series = np.array([wavelength.value ** n for n in range(len(a))])

    vl = np.sum(a * wavelength_series)

    return vl


def _get_limb_dark(sun_map):
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
    DOI : https://doi.org/10.1063/1.1325201
    """

    if not (4000 * u.AA < sun_map.wavelength < 15000 * u.AA):
        raise ValueError("The wavelength of the observation " +
                          "must be between 4000 and 15000 Angstrom.")

    ul = _darklimb_u(sun_map.wavelength)
    vl = _darklimb_v(sun_map.wavelength)

    x_center, y_center = sun_map.world_to_pixel(SkyCoord(0, 0, unit=u.arcsec, frame=sun_map.coordinate_frame))
    x_dim = sun_map.dimensions.x
    y_dim = sun_map.dimensions.y
    radius = sun_map.rsun_obs / sun_map.scale.axis1

    x_grid, y_grid = np.meshgrid(np.arange(x_dim.to_value()), np.arange(y_dim.to_value()), sparse=True)

    x_2 = (x_grid - x_center.value) ** 2
    y_2 = (y_grid - y_center.value) ** 2

    dist_grid = np.sqrt(x_2 + y_2) / radius.value
    dist_grid[np.where(dist_grid > 1)] = 0

    limbfilt = 1 - ul - vl + ul * np.cos(np.arcsin(dist_grid)) + vl * np.cos(np.arcsin(dist_grid)) ** 2

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
    """
    new_data = sun_map.data / _get_limb_dark(sun_map)
    new_meta = copy.deepcopy(sun_map.meta)
    new_meta['history'] += "Limb Darkening effects removed. The limb darkening function \
                            uses a 5th order polynomial fitting to the limb darkening constants obtained from Astrophysical Quantities."

    return (sun_map._new_instance(new_data, new_meta, sun_map.plot_settings))
