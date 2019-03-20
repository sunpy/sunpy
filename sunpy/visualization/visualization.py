"""
This module provides plotting support in iPython.
"""
from matplotlib import pyplot as plt

__all__ = ["toggle_pylab", "axis_labels_from_ctype"]


def toggle_pylab(fn):
    """
    A decorator to prevent functions from opening Matplotlib windows
    unexpectedly when SunPy is run in interactive shells like iPython.

    Toggles the value of `matplotlib.pyplot.isinteractive`.
    """
    if plt.isinteractive():
        def fn_itoggle(*args, **kwargs):
            plt.ioff()
            ret = fn(*args, **kwargs)
            plt.ion()
            return ret
        return fn_itoggle
    else:
        return fn


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
