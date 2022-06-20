"""
This module provides functions that draw on Map plots.
"""
import numpy as np
from matplotlib.lines import Line2D

import astropy.units as u
from astropy.coordinates import SkyCoord

from sunpy.coordinates import HeliographicStonyhurst
from sunpy.visualization import wcsaxes_compat

__all__ = ["draw_equator", "draw_prime_meridian"]


def draw_equator(axes, resolution=300, **kwargs):
    """
    Draws the equator as seen by the axes observer.

    Parameters
    ----------
    axes : `matplotlib.axes.Axes`
        The axes to plot the equator on.
    resolution : `int`
        The number of points used to represent the equator.

    Examples
    --------
    >>> import sunpy.map
    >>> from sunpy.visualization import draw
    >>> import sunpy.data.sample   # doctest: +REMOTE_DATA
    >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)   # doctest: +REMOTE_DATA
    >>> aia.plot()   # doctest: +SKIP
    >>> draw.draw_equator(aia._check_axes(None))   # doctest: +SKIP
    """
    if not wcsaxes_compat.is_wcsaxes(axes):
        raise ValueError('axes must be a WCSAxes')
    axes_frame = axes._transform_pixel2world.frame_out
    observer = axes_frame.observer

    c_kw = {'linewidth': 1,
            'color': 'white',
            'zorder': 100}
    c_kw.update(kwargs)

    transform = axes.get_transform(axes_frame)
    c_kw.setdefault('transform', transform)

    lat = 0*u.deg
    lat0 = SkyCoord(np.linspace(-90, 90, resolution) * u.deg,
                    np.ones(resolution) * lat,
                    frame=HeliographicStonyhurst)
    lat0 = lat0.transform_to(axes_frame)

    rsun = axes_frame.rsun
    reference_dist = np.sqrt(observer.radius**2 - rsun**2)
    is_visible = lat0.spherical.distance <= reference_dist
    vis = np.where(is_visible == True)
    visible = lat0[vis]

    equator = Line2D(visible.Tx.deg, visible.Ty.deg,
                     **c_kw)
    axes.add_line(equator)
    return equator

def draw_prime_meridian(axes, resolution=300, **kwargs):
    """
    Draws the prime meridian as seen by the axes observer.
    Any part that is hidden is not plotted.

    Parameters
    ----------
    axes : `matplotlib.axes.Axes`
        The axes to plot the prime meridian on.
    resolution : `int`
        The number of points used to represent the prime meridian.
    """
    if not wcsaxes_compat.is_wcsaxes(axes):
        raise ValueError('axes must be a WCSAxes')
    axes_frame = axes._transform_pixel2world.frame_out
    observer = axes_frame.observer

    c_kw = {'linewidth': 1,
            'color': 'white',
            'zorder': 100}
    c_kw.update(kwargs)

    transform = axes.get_transform(axes_frame)
    c_kw.setdefault('transform', transform)

    lon = 0*u.deg
    lon0 = SkyCoord(np.ones(resolution) * lon,
                    np.linspace(-90, 90, resolution) * u.deg,
                    frame=HeliographicStonyhurst)
    lon0 = lon0.transform_to(axes_frame)

    rsun = axes_frame.rsun
    reference_dist = np.sqrt(observer.radius**2 - rsun**2)
    is_visible = lon0.spherical.distance <= reference_dist
    vis = np.where(is_visible == True)
    visible = lon0[vis]

    prime_meridian = Line2D(visible.Tx.deg, visible.Ty.deg,
                            **c_kw)
    axes.add_line(prime_meridian)
    return prime_meridian
