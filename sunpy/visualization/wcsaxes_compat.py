# -*- coding: utf-8 -*-
"""
Helpers and Functions to make WCSAxes work in SunPy
"""
import warnings

import matplotlib.pyplot as plt

import astropy.units as u

try:
    import wcsaxes
    HAVE_WCSAXES = True

except ImportError:
    HAVE_WCSAXES = False
    warnings.warn("SunPy plotting is improved by installing the WCSAxes module: http://wcsaxes.readthedocs.io/en/latest/index.html")

FORCE_NO_WCSAXES = False

__all__ = ['HAVE_WCSAXES', 'is_wcsaxes', 'FORCE_NO_WCSAXES']

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

    if HAVE_WCSAXES and not FORCE_NO_WCSAXES:
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
        if HAVE_WCSAXES and not FORCE_NO_WCSAXES:
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

def default_wcs_grid(axes):
    """
    Apply some default wcsaxes grid formatting.

    Parameters
    ----------
    axes : `wcsaxes.WCSAxes` object.
        The `~wcsaxes.WCSAxes` object to draw the world coordinate grid on.
    """
    if not isinstance(axes, wcsaxes.WCSAxes):
        raise TypeError("This axes is not a WCSAxes")

    x = axes.coords[0]
    y = axes.coords[1]

    x.set_ticks(color='white')
    y.set_ticks(color='white')

    x.set_ticks_position('bl')
    y.set_ticks_position('bl')

    if x.coord_type != 'longitude':
        x.set_coord_type('longitude', coord_wrap=180.)
    if y.coord_type != 'latitude':
        y.set_coord_type('latitude')

    x.set_major_formatter('s.s')
    y.set_major_formatter('s.s')

    axes.coords.grid(color='white', alpha=0.6)

def wcsaxes_heliographic_overlay(axes):
    """
    Create a heliographic overlay using wcsaxes.

    Also draw a grid and label the top axes.

    Parameters
    ----------
    axes : `wcsaxes.WCSAxes` object.
        The `~wcsaxes.WCSAxes` object to create the HGS overlay on.

    Returns
    -------
    overlay : wcsaxes overlay
        The overlay object.

    .. note::
        This function requires the WIP sunpy.coordinates module.

    """
    overlay = axes.get_coords_overlay('heliographicstonyhurst')

    lon = overlay[0]
    lat = overlay[1]

    lon.coord_wrap = 180
    lon.set_major_formatter('dd')

    lon.set_axislabel('Solar Longitude')
    lat.set_axislabel('Solar Latitude')

    lon.set_ticks_position('tr')
    lat.set_ticks_position('tr')

    lon.set_ticks(spacing=10. * u.deg, color='white')
    lat.set_ticks(spacing=10. * u.deg, color='white')

    overlay.grid(color='white', alpha=0.5)

    return overlay
