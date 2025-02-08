"""
This module provides plotting support in iPython.
"""
from functools import wraps

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.visualization.wcsaxes import CoordinateHelper

__all__ = ['peek_show', "axis_labels_from_ctype", "show_hpr_impact_angle"]


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
    ctype : `str`
        Coordinate type.
    unit : `str`, `None`
        Required unit. If `None` no unit is added to the label.

    Returns
    -------
    `str`
        "Axis Label [Unit]"
    """
    ctype_short = ctype[:4]

    labels = {'HGLN': 'Heliographic Longitude',
              'CRLN': 'Carrington Longitude',
              'HPLN': 'Helioprojective Longitude (Solar-X)',
              'HRLN': 'Helioprojective Position Angle',
              'SOLX': 'Heliocentric X',

              'HGLT': 'Latitude',
              'CRLT': 'Latitude',
              'HPLT': 'Helioprojective Latitude (Solar-Y)',
              'HRLT': 'Helioprojective Declination',
              'SOLY': 'Heliocentric Y'}

    label = labels.get(ctype_short, f"{ctype}")
    if unit is not None:
        label += f' [{unit}]'

    return label


def show_hpr_impact_angle(declination_axis):
    """
    Modify a declination axis for Helioprojective Radial plots to show impact angle
    instead.

    In order to accommodate the FITS-WCS machinery, the declination is used instead
    of the impact angle (which is the declination plus 90 degrees). This function
    changes the tick-label formatter for an axis so that it adds 90 degrees before
    rendering the labels.

    Parameters
    ----------
    coordinate_axis : `astropy.visualization.wcsaxes.CoordinateHelper`
        The coordinate axis for Helioprojective Radial declination (``delta``)

    See Also
    --------
    sunpy.coordinates.HelioprojectiveRadial

    Examples
    --------
    .. minigallery:: sunpy.visualization.show_hpr_impact_angle
    """
    if not isinstance(declination_axis, CoordinateHelper):
        raise TypeError("The input should be one of the WCSAxes coordinate axes.")

    old_format_method = declination_axis.formatter

    def new_format_method(values, spacing, format="auto"):
        return old_format_method(values + 90*u.deg, spacing, format=format)
    # TODO: Astropy >= 7.0 use declination_axis.set_major_formatter(new_format_method)
    declination_axis._formatter_locator.formatter = new_format_method
    declination_axis.set_axislabel('Helioprojective Impact Angle')
