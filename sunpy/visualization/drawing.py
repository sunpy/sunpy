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
from astropy.wcs.wcsapi import BaseHighLevelWCS

from sunpy.coordinates import HeliographicCarrington, HeliographicStonyhurst
from sunpy.coordinates.frames import HeliocentricInertial, Helioprojective
from sunpy.coordinates.sun import _angular_radius
from sunpy.coordinates.utils import get_limb_coordinates
from sunpy.util import grid_perimeter
from sunpy.visualization import wcsaxes_compat

__all__ = ["limb", "equator", "prime_meridian", "extent"]


@u.quantity_input
def limb(axes, observer, *, rsun: u.m = R_sun, resolution=1000, **kwargs):
    """
    Draws the solar limb as seen by the specified observer.

    The limb is a circle for only the simplest plots. If the specified
    observer of the limb is different from the observer of the coordinate frame
    of the plot axes, not only may the limb not be a true circle, a portion of
    the limb may be hidden from the observer. In that case, the circle is
    divided into visible and hidden segments, represented by solid and dotted
    lines, respectively.

    Parameters
    ----------
    axes : `~matplotlib.axes` or ``None``
        Axes to plot limb on.
    observer : `astropy.coordinates.SkyCoord`
        Observer coordinate for which the limb is drawn.
    rsun : `~astropy.units.Quantity`
        Solar radius (in physical length units) at which to draw the limb.
        Defaults to the standard photospheric radius.
    resolution : `int`
        The number of points to use to represent the limb.

    Returns
    -------
    visible : `~matplotlib.patches.Polygon` or `~matplotlib.patches.Circle` or None
        The patch added to the axes for the visible part of the limb (i.e., the
        "near" side of the Sun).
    hidden : `~matplotlib.patches.Polygon` or None
        The patch added to the axes for the hidden part of the limb (i.e., the
        "far" side of the Sun).

    Notes
    -----
    Keyword arguments are passed onto the patches.

    If the limb is a true circle, ``visible`` will instead be
    `~matplotlib.patches.Circle` and ``hidden`` will be ``None``.

    If there are no hidden points (e.g., on a synoptic map any limb is fully
    visible) ``hidden`` will be ``None``.

    If there are no visible points (e.g., for an observer on the opposite side
    of the Sun to the map observer) ``visible`` will be ``None``.

    To avoid triggering Matplotlib auto-scaling, these patches are added as
    artists instead of patches. One consequence is that the plot legend is not
    populated automatically when the limb is specified with a text label. See
    :ref:`sphx_glr_gallery_text_labels_and_annotations_custom_legends.py` in
    the Matplotlib documentation for examples of creating a custom legend.
    """
    if not wcsaxes_compat.is_wcsaxes(axes):
        raise ValueError('axes must be a WCSAxes')

    c_kw = {'fill': False,
            'color': 'white',
            'zorder': 100}
    c_kw.update(kwargs)

    transform = axes.get_transform('world')
    # transform is always passed on as a keyword argument
    c_kw.setdefault('transform', transform)

    # If the observer matches the axes's observer frame and is Helioprojective, use a Circle
    axes_frame = wcsapi_to_celestial_frame(axes.wcs)
    if isinstance(axes_frame, Helioprojective):
        axes_observer = SkyCoord(axes_frame.observer)
        if axes_observer.separation_3d(observer) < 1 * u.m:
            distance = observer.transform_to(HeliocentricInertial).spherical.distance
            # Obtain the solar radius and the world->pixel transform
            angular_radius = _angular_radius(rsun, distance)
            circ = patches.Circle([0, 0], radius=angular_radius.to_value(u.deg), **c_kw)
            axes.add_artist(circ)
            return circ, None

    # Otherwise, we use Polygon to be able to distort the limb
    # Create the limb coordinate array using Heliocentric Radial
    limb = get_limb_coordinates(observer, rsun, resolution)

    # Transform the limb to the axes frame and get the 2D vertices
    visible, hidden = _plot_vertices(limb, axes, axes_frame, rsun, **kwargs)

    return visible, hidden


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
    >>> from sunpy.visualization import drawing
    >>> import sunpy.data.sample   # doctest: +REMOTE_DATA
    >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)   # doctest: +REMOTE_DATA +IGNORE_WARNINGS
    >>> fig = plt.figure()   # doctest: +SKIP
    >>> ax = fig.add_subplot(projection=aia)   # doctest: +SKIP
    >>> aia.plot()   # doctest: +SKIP
    >>> drawing.equator(ax)   # doctest: +SKIP
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


@u.quantity_input
def extent(axes, wcs, **kwargs):
    """
    Draws the extent as defined by a given `~astropy.wcs.WCS` instance.

    Parameters
    ----------
    axes : `astropy.visualization.wcsaxes.WCSAxes`
        The axes to plot the WCS extent on, or "None" to use current axes.
    wcs : `~astropy.wcs.wcsapi.BaseLowLevelWCS`
        The WCS that defines the extent to be drawn.

    Returns
    -------
    visible : `~matplotlib.patches.Polygon`
        The patch added to the axes for the visible part of the WCS extent.
    hidden : `~matplotlib.patches.Polygon`
        The patch added to the axes for the hidden part of the WCS extent.
    """
    if not wcsaxes_compat.is_wcsaxes(axes):
        raise ValueError('axes must be a WCSAxes')
    if not isinstance(wcs, BaseHighLevelWCS):
        raise TypeError("wcs should be a High Level WCS object")
    axes_frame = wcsapi_to_celestial_frame(axes.wcs)
    # Traverse the edges of the WCS in pixel space
    edge_pixels = grid_perimeter(*wcs.low_level_wcs.pixel_shape) - 0.5
    edge_coords = wcs.pixel_to_world(*edge_pixels.T)
    # Filter out any non-SkyCoord coordinates returned for these pixel axes
    if not isinstance(edge_coords, SkyCoord):
        edge_coords = [c for c in edge_coords if isinstance(c, SkyCoord)][0]
    visible, hidden = _plot_vertices(edge_coords, axes, axes_frame, R_sun, close_path=True, **kwargs)
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

    # Determine which points are visible (2D points are always visible)
    is_2d = (norm := coord.spherical.norm()).unit is u.one and u.allclose(norm, 1*u.one)
    if not is_2d and hasattr(frame, 'observer'):
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
    if close_path is not True:
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
