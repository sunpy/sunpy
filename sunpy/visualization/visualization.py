"""
This module provides plotting support in iPython.
"""
from functools import wraps

from matplotlib import pyplot as plt

__all__ = ['peek_show', "axis_labels_from_ctype"]


def peek_show(func):
    """
    A decorator to place on ``peek()`` methods to show the figure.

    The ``peek()`` method should return the figure then this method will
    attempt to show it in the correct way. This decorator will not return the
    figure to the user.
    """
    @wraps(func)
    def show_if_interactive(*args, **kwargs):
        figure = func(*args, **kwargs)
        # Show the figure if using an interactive Matplotlib backend
        if mpl.get_backend() in mpl.rcsetup.interactive_bk:
            figure.show()

        # NOTE: We do not return `figure` here because `peek()` methods return `None`.

    return show_if_interactive


def axis_labels_from_ctype(ctype, unit):
    """
    Returns axis labels for the given coordinate type and unit.

    Parameters
    ----------
    ctype: `str`
        Coordinate type.
    unit: `str`
        Required unit.

    Returns
    -------
    `str`
        "Axis Label [Unit]"
    """
    ctype_short = ctype[:4]

    labels = {'HGLN': 'Heliographic Longitude [{}]'.format(unit),
              'CRLN': 'Carrington Longitude [{}]'.format(unit),
              'HPLN': 'Helioprojective Longitude (Solar-X) [{}]'.format(unit),
              'SOLX': 'Heliocentric X [{}]'.format(unit),

              'HGLT': 'Latitude [{}]'.format(unit),
              'CRLT': 'Latitude [{}]'.format(unit),
              'HPLT': 'Helioprojective Latitude (Solar-Y) [{}]'.format(unit),
              'SOLY': 'Heliocentric Y [{}]'.format(unit)}

    return labels.get(ctype_short, "{} [{}]".format(ctype, unit))
