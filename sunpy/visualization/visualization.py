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

    labels = {'HGLN': f'Heliographic Longitude [{unit}]',
              'CRLN': f'Carrington Longitude [{unit}]',
              'HPLN': f'Helioprojective Longitude (Solar-X) [{unit}]',
              'SOLX': f'Heliocentric X [{unit}]',

              'HGLT': f'Latitude [{unit}]',
              'CRLT': f'Latitude [{unit}]',
              'HPLT': f'Helioprojective Latitude (Solar-Y) [{unit}]',
              'SOLY': f'Heliocentric Y [{unit}]'}

    return labels.get(ctype_short, f"{ctype} [{unit}]")
