"""
This module provides functions that draw on Astropy's WCSAxes.
"""
import numpy as np
from matplotlib import patches
from matplotlib.path import Path

import astropy.units as u
from astropy.constants import R_sun
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes.wcsapi import wcsapi_to_celestial_frame

from sunpy.coordinates import HeliographicCarrington, HeliographicStonyhurst
from sunpy.visualization import wcsaxes_compat

__all__ = ["equator", "prime_meridian"]


@u.quantity_input
def equator(axes, *, rsun: u.m = R_sun, resolution=500, **kwargs):
    """
    Draws the solar equator as seen by the axes observer.
    Hidden parts are drawn as a dotted line.

    Parameters
    ----------
    axes : `matplotlib.axes.Axes`
        The axes to plot the equator on.
    rsun : `~astropy.units.Quantity`
        Solar radius (in physical length units) at which to draw the solar
        equator. Defaults to the standard photospheric radius.
    resolution : `int`
        The number of points used to represent the equator.

    Returns
    -------
    visible : `~matplotlib.patches.Polygon`
        The patch added to the axes for the visible part of the solar equator.
    hidden : `~matplotlib.patches.Polygon`
        The patch added to the axes for the hidden part of the solar equator.

    Examples
    --------
    >>> import sunpy.map
    >>> from sunpy.visualization import draw
    >>> import sunpy.data.sample   # doctest: +REMOTE_DATA
    >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)   # doctest: +REMOTE_DATA
    >>> fig = plt.figure()   # doctest: +SKIP
    >>> ax = fig.add_subplot(projection=aia)   # doctest: +SKIP
    >>> aia.plot()   # doctest: +SKIP
    >>> draw.equator(ax)   # doctest: +SKIP
    """
    if not wcsaxes_compat.is_wcsaxes(axes):
        raise ValueError('axes must be a WCSAxes')
    axes_frame = wcsapi_to_celestial_frame(axes.wcs)

    lat = 0*u.deg
    lat0 = SkyCoord(np.linspace(-180, 179, resolution) * u.deg,
                    np.ones(resolution) * lat, radius=rsun,
                    frame=HeliographicStonyhurst(obstime=axes_frame.obstime))
    visible, hidden = _plot_vertices(lat0, axes, axes_frame, rsun,
                                     **kwargs)
    return visible, hidden


@u.quantity_input
def prime_meridian(axes, *, rsun: u.m = R_sun, resolution=500, **kwargs):
    """
    Draws the solar prime meridian (zero Carrington longitude) as seen by the
    axes observer. Hidden parts are drawn as a dotted line.

    Parameters
    ----------
    axes : `matplotlib.axes.Axes`
        The axes to plot the prime meridian on, or "None" to use current axes.
    rsun : `~astropy.units.Quantity`
        Solar radius (in physical length units) at which to draw the solar
        prime meridian. Defaults to the standard photospheric radius.
    resolution : `int`
        The number of points used to represent the prime meridian.

    Returns
    -------
    visible : `~matplotlib.patches.Polygon`
        The patch added to the axes for the visible part of the solar equator.
    hidden : `~matplotlib.patches.Polygon`
        The patch added to the axes for the hidden part of the solar equator.
    """
    if not wcsaxes_compat.is_wcsaxes(axes):
        raise ValueError('axes must be a WCSAxes')
    axes_frame = wcsapi_to_celestial_frame(axes.wcs)

    if not hasattr(axes_frame, 'observer'):
        raise ValueError('the coordinate frame of the WCSAxes does not have an observer, '
                         'so zero Carrington longitude cannot be determined.')
    observer = axes_frame.observer

    lon = 0*u.deg
    lon0 = SkyCoord(np.ones(resolution) * lon,
                    np.linspace(-90, 90, resolution) * u.deg, radius=rsun,
                    frame=HeliographicCarrington(observer=observer,
                                                 obstime=axes_frame.obstime))
    visible, hidden = _plot_vertices(lon0, axes, axes_frame, rsun,
                                     close_path=False, **kwargs)
    return visible, hidden


def _plot_vertices(coord, axes, frame, rsun, close_path=True, **kwargs):
    """
    Draws the provided SkyCoord on the WCSAxes as `~matplotlib.patches.Polygon`
    objects depending on visibility.
    """
    c_kw = {'fill': False,
            'color': 'white',
            'zorder': 100}
    c_kw.update(kwargs)
    transform = axes.get_transform("world")
    c_kw.setdefault('transform', transform)

    # Get the 2D vertices of the coordinates
    coord = coord.transform_to(frame)
    Tx = coord.spherical.lon.to_value(u.deg)
    Ty = coord.spherical.lat.to_value(u.deg)
    vertices = np.array([Tx, Ty]).T
    # Determine which points are visible
    if hasattr(frame, 'observer'):
        # The reference distance is the distance to the limb for the axes
        # observer
        rsun = getattr(frame, 'rsun', rsun)
        reference_distance = np.sqrt(frame.observer.radius**2 - rsun**2)
        is_visible = coord.spherical.distance <= reference_distance
    else:
        # If the axes has no observer, the entire limb is considered visible
        is_visible = np.ones_like(coord.spherical.distance, bool, subok=False)

    # Identify discontinuities. Uses the same approach as
    # astropy.visualization.wcsaxes.grid_paths.get_lon_lat_path()
    step = np.sqrt((vertices[1:, 0] - vertices[:-1, 0]) ** 2 +
                   (vertices[1:, 1] - vertices[:-1, 1]) ** 2)
    continuous = np.concatenate([[True, True], step[1:] < 100 * step[:-1]])
    if not close_path is True:
        continuous[0] = False
        continuous[-1] = False

    visible, hidden = None, None
    if np.sum(is_visible) > 0:
        if 'linestyle' not in kwargs:
            c_kw['linestyle'] = '-'
        visible = patches.Polygon(vertices, **c_kw)
        _modify_polygon_visibility(visible, is_visible & continuous)
        axes.add_artist(visible)
    if np.sum(~is_visible) > 0:
        if 'linestyle' not in kwargs:
            c_kw['linestyle'] = ':'
        hidden = patches.Polygon(vertices, **c_kw)
        _modify_polygon_visibility(hidden, ~is_visible & continuous)
        axes.add_artist(hidden)

    return visible, hidden


def _modify_polygon_visibility(polygon, keep):
    polygon_codes = polygon.get_path().codes
    polygon_codes[:-1][~keep] = Path.MOVETO
    polygon_codes[-1] = Path.MOVETO if not keep[0] else Path.LINETO
