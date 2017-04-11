# -*- coding: utf-8 -*-
"""
Helpers and Functions to make WCSAxes work in SunPy
"""
import matplotlib.pyplot as plt

import astropy.units as u

from astropy.utils.misc import isiterable

try:
    from astropy.visualization import wcsaxes
except ImportError:
    raise ImportError("Astropy >= 1.3 is required to use SunPy")

#  Force is put here to enable disabling all checks in this module. It should
#  only be used by tests and other such hacks.
_FORCE_NO_WCSAXES = False

__all__ = ['is_wcsaxes']


def is_wcsaxes(axes):
    """
    Test a matplotlib Axes object to see if it is an instance of WCSAxes.

    Parameters
    ----------
    axes : `matplotlib.axes` Object
        Axes to test

    Returns
    -------
    result : `bool`
        Result of the test
    """

    if not _FORCE_NO_WCSAXES:
        return isinstance(axes, wcsaxes.WCSAxes)
    else:
        return False


def gca_wcs(wcs, fig=None):
    """
    Get the current axes, and return a WCSAxes if possible.

    Parameters
    ----------
    wcs : `astropy.wcs.WCS`
        A `~astropy.wcs.WCS` object used to create a new axes.
    fig : `matplotlib.figure.Figure`
        The figure in which to check for the axes.

    Returns
    -------
    ax : `matplotlib.axes.Axes` or `wcsaxes.WCSAxes` object.
        The current axes, or a new one if created.

    """

    if not fig:
        fig = plt.gcf()

    if not len(fig.get_axes()):
        if not _FORCE_NO_WCSAXES:
            ax = plt.gca(projection=wcs)
        else:
            ax = plt.gca()

    else:
        ax = plt.gca()

    return ax


def get_world_transform(axes):
    """
    Get the transformation to world coordinates.

    If the axes is a `wcaxes.WCSAxes` instance this returns the transform to
    the ``'world'`` coordinates, otherwise it returns the transform to the
    matplotlib data coordinates, which are assumed to be in world coordinates.

    Parameters
    ----------
    axes : `wcsaxes.WCSAxes` or `matplotlib.axes.Axes` obejct.
        The axes to get the transform from.

    Returns
    -------
    transform : `matplotlib.transforms.CompositeGenericTransform`
        The transformation object.
    """
    if is_wcsaxes(axes):
        transform = axes.get_transform('world')
    else:
        transform = axes.transData

    return transform


def ctype_longitude_wrap(ctype):
    """
    Return the longitude wrapping for a given ctype. If the wrapping is either
    360 or not specified then this function returns `None`.
    """
    if ctype[:4] in ['HPLN', 'HGLN']:
        return 180.


def solar_coord_type_from_ctype(ctype):
    """
    Determine whether a particular WCS ctype corresponds to an angle or scalar
    coordinate.
    """
    wrapping = ctype_longitude_wrap(ctype)

    if ctype[2:4] == 'LN':
        return 'longitude', wrapping
    elif ctype[2:4] == 'LT':
        return 'latitude', wrapping
    else:
        return 'scalar', None


def default_wcs_grid(axes, units, ctypes):
    """
    Apply some default wcsaxes grid formatting.

    Parameters
    ----------
    axes : `wcsaxes.WCSAxes` object.
        The `~wcsaxes.WCSAxes` object to draw the world coordinate grid on.

    units : `tuple`
        The axes units axes x y order.
    """
    if not isinstance(axes, wcsaxes.WCSAxes):
        raise TypeError("This axes is not a WCSAxes")

    x = axes.coords[0]
    y = axes.coords[1]

    x.set_ticks(color='white')
    y.set_ticks(color='white')

    x.set_ticks_position('bl')
    y.set_ticks_position('bl')

    x.set_coord_type(*solar_coord_type_from_ctype(ctypes[0]))
    y.set_coord_type(*solar_coord_type_from_ctype(ctypes[1]))

    if units[0] is u.deg:
        x.set_major_formatter('dd')
    elif units[0] is u.arcsec:
        x.set_major_formatter('s.s')
    else:
        x.set_major_formatter('x.x')

    if units[1] is u.deg:
        y.set_major_formatter('dd')
    elif units[1] is u.arcsec:
        y.set_major_formatter('s.s')
    else:
        y.set_major_formatter('x.x')

    axes.coords.grid(color='white', alpha=0.6, linestyle='dotted',
                     linewidth=0.5)


@u.quantity_input(grid_spacing=u.deg)
def wcsaxes_heliographic_overlay(axes, grid_spacing=10*u.deg):
    """
    Create a heliographic overlay using wcsaxes.

    Also draw a grid and label the top axes.

    Parameters
    ----------
    axes : `wcsaxes.WCSAxes` object.
        The `~wcsaxes.WCSAxes` object to create the HGS overlay on.

    grid_spacing: `astropy.units.Quantity`
        Spacing for longitude and latitude grid in degrees.

    Returns
    -------
    overlay : wcsaxes overlay
        The overlay object.
    """

    # Unpack spacing
    if isinstance(grid_spacing, u.Quantity) and grid_spacing.size == 1:
        lon_space = lat_space = grid_spacing
    elif grid_spacing.size == 2:
        lon_space, lat_space = grid_spacing
    else:
        raise ValueError("grid_spacing must be a Quantity of length one or two.")

    overlay = axes.get_coords_overlay('heliographic_stonyhurst')

    lon = overlay[0]
    lat = overlay[1]

    lon.coord_wrap = 180
    lon.set_major_formatter('dd')

    lon.set_axislabel('Solar Longitude')
    lat.set_axislabel('Solar Latitude')

    lon.set_ticks_position('tr')
    lat.set_ticks_position('tr')

    lon.set_ticks(spacing=lon_space, color='white')
    lat.set_ticks(spacing=lat_space, color='white')

    overlay.grid(color='white', alpha=0.5)

    if axes.title:
        x, y = axes.title.get_position()
        axes.title.set_position([x, y + 0.05])

    return overlay
