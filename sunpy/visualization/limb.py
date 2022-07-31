import astropy.units as u
from astropy.constants import R_sun

from sunpy.util.decorators import deprecated
from sunpy.visualization import drawing

__all__ = ['draw_limb']


@deprecated(since="4.1", message="This will be moved.",
            alternative="sunpy.visualization.drawing.limb")
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
    """
    visible, hidden = drawing.limb(axes, observer, rsun=rsun,
                                   resolution=resolution, **kwargs)
    return visible, hidden
