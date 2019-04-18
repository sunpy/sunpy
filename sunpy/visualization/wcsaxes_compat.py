"""
This module provides functions to make WCSAxes work in SunPy.
"""
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.visualization import wcsaxes

# Force is put here to enable disabling all checks in this module.
# It should only be used by tests and other such hacks.
_FORCE_NO_WCSAXES = False

__all__ = ["is_wcsaxes", "gca_wcs", "get_world_transform",
           "default_wcs_grid", "wcsaxes_heliographic_overlay"]


def is_wcsaxes(axes):
    """
    Tests a `matplotlib.axes.Axes` object to see if it is an instance of
    `~astropy.visualization.wcsaxes.WCSAxes`.

    Parameters
    ----------
    axes : `matplotlib.axes`
        Axes to test.

    Returns
    -------
    `bool`
        Result of the test.
    """
    if not _FORCE_NO_WCSAXES:
        return isinstance(axes, wcsaxes.WCSAxes)
    else:
        return False


def gca_wcs(wcs, fig=None, slices=None):
    """
    Get the current axes, and return a `~astropy.visualization.wcsaxes.WCSAxes`
    if possible.

    Parameters
    ----------
    wcs : `astropy.wcs.WCS`
        A `~astropy.wcs.WCS` object used to create a new axes.
    fig : `matplotlib.figure.Figure`
        The figure in which to check for the axes.
    slices : `tuple`
        ``slices`` is passed to `~astropy.visualization.wcsaxes.WCSAxes` to describe
        which two dimensions of the `~astropy.wcs.WCS` object are being plotted.
        This slices the multidimensional wcs object in the way it needs to be sliced.

    Returns
    -------
    `matplotlib.axes.Axes` or `~astropy.visualization.wcsaxes.WCSAxes`
        The current axes, or a new one if created.
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
    returns the transform to the "world" coordinates, otherwise it returns
    the transform to the matplotlib data coordinates, which are assumed to be in
    world coordinates.

    Parameters
    ----------
    axes : `~astropy.visualization.wcsaxes.WCSAxes` or `~matplotlib.axes.Axes`
        The axes to get the transform from.

    Returns
    -------
    `~matplotlib.transforms.CompositeGenericTransform`
        The transformation object.
    """
    if is_wcsaxes(axes):
        transform = axes.get_transform('world')
    else:
        transform = axes.transData

    return transform


def default_wcs_grid(axes):
    """
    Apply some default `~astropy.visualization.wcsaxes.WCSAxes` grid
    formatting.

    Parameters
    ----------
    axes : `~astropy.visualization.wcsaxes.WCSAxes`
        The `~astropy.visualization.wcsaxes.WCSAxes` object to draw the world
        coordinate grid on.
    """
    axes.coords.grid(color='white', alpha=0.6, linestyle='dotted',
                     linewidth=0.5)


@u.quantity_input
def wcsaxes_heliographic_overlay(axes, grid_spacing: u.deg = 10*u.deg, **kwargs):
    """
    Create a heliographic overlay using
    `~astropy.visualization.wcsaxes.WCSAxes`.

    Will draw a grid and label the top axes.

    Parameters
    ----------
    axes : `~astropy.visualization.wcsaxes.WCSAxes`
        The `~astropy.visualization.wcsaxes.WCSAxes` object to create the HGS overlay on.
    grid_spacing: `~astropy.units.Quantity`
        Spacing for longitude and latitude grid in degrees.

    Returns
    -------
    `~astropy.visualization.wcsaxes.WCSAxes`
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
