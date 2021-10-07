"""
This module provides plotting support in iPython.
"""
import sunpy.util.visualization as viz
from sunpy.util.decorators import deprecated

__all__ = ['peek_show', "axis_labels_from_ctype"]


@deprecated('3.1')
def peek_show(func):
    """
    A decorator to place on ``peek()`` methods to show the figure.

    The ``peek()`` method should return the figure then this method will
    attempt to show it in the correct way. This decorator will not return the
    figure to the user.
    """
    return viz.peek_show(func)


@deprecated('3.1')
def axis_labels_from_ctype(ctype, unit):
    """
    Returns axis labels for the given coordinate type and unit.

    Parameters
    ----------
    ctype: `str`
        Coordinate type.
    unit: `str`, `None`
        Required unit. If `None` no unit is added to the label.

    Returns
    -------
    `str`
        "Axis Label [Unit]"
    """
    return viz.axis_labels_from_ctype(ctype, unit)
