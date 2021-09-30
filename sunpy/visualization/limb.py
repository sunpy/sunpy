import numpy as np
from matplotlib import patches
from matplotlib.path import Path

import astropy.units as u
from astropy.constants import R_sun
from astropy.coordinates import SkyCoord

from sunpy.coordinates.frames import HeliocentricInertial, Helioprojective
from sunpy.coordinates.sun import _angular_radius
from sunpy.coordinates.utils import get_limb_coordinates
from sunpy.visualization import wcsaxes_compat

__all__ = ['draw_limb']


@u.quantity_input
def draw_limb(axes, observer, *, rsun: u.m = R_sun, resolution=1000, **kwargs):
    """
    Draws the solar limb as seen by the specified observer.

    The limb is a circle for only the simplest plots.  If the specified
    observer of the limb is different from the observer of the coordinate frame
    of the plot axes, not only may the limb not be a true circle, a portion of
    the limb may be hidden from the observer.  In that case, the circle is
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
    artists instead of patches.  One consequence is that the plot legend is not
    populated automatically when the limb is specified with a text label.  See
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
    axes_frame = axes._transform_pixel2world.frame_out
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
    limb_in_axes = limb.transform_to(axes_frame)
    Tx = limb_in_axes.spherical.lon.to_value(u.deg)
    Ty = limb_in_axes.spherical.lat.to_value(u.deg)
    vertices = np.array([Tx, Ty]).T

    # Determine which points are visible
    if hasattr(axes_frame, 'observer'):
        # The reference distance is the distance to the limb for the axes
        # observer
        rsun = getattr(axes_frame, 'rsun', rsun)
        reference_distance = np.sqrt(axes_frame.observer.radius**2 - rsun**2)
        is_visible = limb_in_axes.spherical.distance <= reference_distance
    else:
        # If the axes has no observer, the entire limb is considered visible
        is_visible = np.ones_like(
            limb_in_axes.spherical.distance, bool, subok=False)

    # Identify discontinuities in the limb. Uses the same approach as
    # astropy.visualization.wcsaxes.grid_paths.get_lon_lat_path()
    step = np.sqrt((vertices[1:, 0] - vertices[:-1, 0]) ** 2 +
                   (vertices[1:, 1] - vertices[:-1, 1]) ** 2)
    continuous = np.concatenate([[True, True], step[1:] < 100 * step[:-1]])

    visible, hidden = None, None
    if np.sum(is_visible) > 0:
        # Create the Polygon for the near side of the Sun (using a solid line)
        if 'linestyle' not in kwargs:
            c_kw['linestyle'] = '-'
        visible = patches.Polygon(vertices, **c_kw)
        _modify_polygon_visibility(visible, is_visible & continuous)
        # Add patches as artists rather than patches to avoid triggering auto-scaling
        axes.add_artist(visible)

    if np.sum(~is_visible) > 0:
        # Create the Polygon for the far side of the Sun (using a dotted line)
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
