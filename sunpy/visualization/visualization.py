"""
This module provides plotting support in iPython.
"""
from functools import wraps

import matplotlib.pyplot as plt

__all__ = ['peek_show', "axis_labels_from_ctype"]


def peek_show(func):
    """
    A decorator to place on ``peek()`` methods to show the figure.

    The ``peek()`` method should return the figure then this method will
    attempt to show it in the correct way. This decorator will not return the
    figure to the user.
    """
    @wraps(func)
    def show_figure(*args, **kwargs):
        _ = func(*args, **kwargs)
        plt.show()

    return show_figure


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
