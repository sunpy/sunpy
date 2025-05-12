"""
This module provides plotting support in iPython.
"""
from functools import wraps

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import Transform

import astropy.units as u
from astropy.visualization.wcsaxes import CoordinateHelper
from astropy.wcs.utils import pixel_to_pixel

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


# TODO: Consider making this transform class public
class _PrecomputedPixelCornersTransform(Transform):
    """
    A matplotlib transform that precomputes the pixel->pixel transformation from a
    given data WCS to the WCS of the given WCSAxes, and stores it as a lookup table
    to avoid redundant coordinate transformations. This transformation is computed
    only for the pixel corners of the data WCS (e.g., (-0.5, -0.5) for the bottom-
    left corner), so the transform may not return desired values when called with
    pixel coordinates other than at the corners of pixels.
    """
    input_dims = 2
    output_dims = 2

    def __init__(self, axes, data_wcs):
        super().__init__()
        self.data_y, self.data_x = np.indices(np.array(data_wcs.array_shape) + 1) - 0.5
        self.axes_x, self.axes_y = pixel_to_pixel(data_wcs, axes.wcs, self.data_x, self.data_y)

    def transform_non_affine(self, values):
        ix, iy = np.ceil(values).astype(int).T
        return np.stack([self.axes_x[iy, ix], self.axes_y[iy, ix]], axis=1)
