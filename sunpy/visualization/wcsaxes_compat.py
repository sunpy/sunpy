# -*- coding: utf-8 -*-
"""
Helpers and Functions to make WCSAxes work in SunPy
"""
import matplotlib.pyplot as plt

import astropy.units as u

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


def gca_wcs(wcs, fig=None, slices=None):
    """
    Get the current axes, and return a WCSAxes if possible.

    Parameters
    ----------
    wcs : `astropy.wcs.WCS`
        A `~astropy.wcs.WCS` object used to create a new axes.
    fig : `matplotlib.figure.Figure`
        The figure in which to check for the axes.
    slices : `tuple`
        ``slices`` is passed to `~astropy.visualization.wcsaxes.WCSAxes`
        to describe which two dimensions of the `~astropy.wcs.WCS` object
        are being plotted.
        This slices the multidimensional wcs object in the way it needs
        to be sliced.

    Returns
    -------
    ax : `matplotlib.axes.Axes` or `~astropy.visualization.wcsaxes.WCSAxes`
        object. The current axes, or a new one if created.

    """

    if not fig:
        fig = plt.gcf()

    if not len(fig.get_axes()):
        if not _FORCE_NO_WCSAXES:
            ax = plt.gca(projection=wcs, slices=slices)
        else:
            ax = plt.gca()

    else:
        ax = plt.gca()

    return ax


def get_world_transform(axes):
    """
    Get the transformation to world coordinates.

    If the axes is a `~astropy.visualization.wcsaxes.WCSAxes` instance this
    returns the transform to the ``'world'`` coordinates, otherwise it returns
    the transform to the matplotlib data coordinates, which are assumed to be in
    world coordinates.

    Parameters
    ----------
    axes : `~astropy.visualization.wcsaxes.WCSAxes` or `~matplotlib.axes.Axes`
        object. The axes to get the transform from.

    Returns
    -------
    transform : `~matplotlib.transforms.CompositeGenericTransform`
        The transformation object.
    """
    if is_wcsaxes(axes):
        transform = axes.get_transform('world')
    else:
        transform = axes.transData

    return transform


def solar_coord_type_from_ctype(ctype):
    """
    Determine whether a particular WCS ctype corresponds to an angle or scalar
    coordinate.
    """

    if ctype[2:4] == 'LN':
        if ctype[:4] in ['HPLN', 'HGLN']:
            return 'longitude', 180.

        return 'longitude', None

    elif ctype[2:4] == 'LT':
        return 'latitude', None

    else:
        return 'scalar', None


def default_wcs_ticks(axes, units, ctypes):
    """
    Set the ticks and axes type on a solar WCSAxes plot.
    """

    if not isinstance(axes, wcsaxes.WCSAxes):
        raise TypeError("This axes is not a WCSAxes")

    x = axes.coords[0]
    y = axes.coords[1]

    if x.ticks.get_tick_out() == 'in':
        x.set_ticks(color='white')
    if y.ticks.get_tick_out() == 'in':
        y.set_ticks(color='white')

    x.set_ticks_position('bl')
    y.set_ticks_position('bl')

    xtype = solar_coord_type_from_ctype(ctypes[0])
    ytype = solar_coord_type_from_ctype(ctypes[1])

    x.set_coord_type(*xtype)
    y.set_coord_type(*ytype)

    if xtype[0] == 'scalar':
        x.set_major_formatter('x.x')
    elif units[0] is u.deg:
        x.set_major_formatter('d.d')
    elif units[0] is u.arcsec:
        x.set_major_formatter('s.s')
    else:
        x.set_major_formatter('x.x')

    if ytype[0] == 'scalar':
        x.set_major_formatter('x.x')
    elif units[1] is u.deg:
        y.set_major_formatter('d.d')
    elif units[1] is u.arcsec:
        y.set_major_formatter('s.s')
    else:
        y.set_major_formatter('x.x')


def default_wcs_grid(axes, units, ctypes):
    """
    Apply some default wcsaxes grid formatting.

    Parameters
    ----------
    axes : `~astropy.visualization.wcsaxes.WCSAxes` object.
        The `~astropy.visualization.wcsaxes.WCSAxes` object to draw the world
        coordinate grid on.

    units : `tuple`
        The axes units axes x y order.
    """

    default_wcs_ticks(axes, units, ctypes)

    axes.coords.grid(color='white', alpha=0.6, linestyle='dotted',
                     linewidth=0.5)


@u.quantity_input(grid_spacing=u.deg)
def wcsaxes_heliographic_overlay(axes, grid_spacing=10*u.deg, **kwargs):
    """
    Create a heliographic overlay using wcsaxes.

    Also draw a grid and label the top axes.

    Parameters
    ----------
    axes : `~astropy.visualization.wcsaxes.WCSAxes` object.
        The `~astropy.visualization.wcsaxes.WCSAxes` object to create the HGS overlay on.

    grid_spacing: `~astropy.units.Quantity`
        Spacing for longitude and latitude grid in degrees.

    Returns
    -------
    overlay : `~astropy.visualization.wcsaxes.WCSAxes` overlay
        The overlay object.

    Notes
    -----
    Keywords are passed to `~astropy.visualization.wcsaxes.coordinates_map.CoordinatesMap.grid`.

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

    lon.set_axislabel('Solar Longitude', minpad=0.8)
    lat.set_axislabel('Solar Latitude', minpad=0.9)

    lon.set_ticks_position('tr')
    lat.set_ticks_position('tr')

    grid_kw = {'color': 'white', 'zorder': 100, 'alpha': 0.5}
    grid_kw.update(kwargs)

    lon.set_ticks(spacing=lon_space, color=grid_kw['color'])
    lat.set_ticks(spacing=lat_space, color=grid_kw['color'])

    overlay.grid(**grid_kw)

    if axes.title:
        x, y = axes.title.get_position()
        axes.title.set_position([x, y + 0.08])

    return overlay
