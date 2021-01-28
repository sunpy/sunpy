# Vendored from astropy.visualization.wcsaxes.patches
# Licensed under a 3-clause BSD style license - see Astropy's LICENSE.rst
#
# This file can be removed when our minimum Astropy dependency is >= 4.2

import numpy as np
from matplotlib.patches import Polygon

from astropy import units as u

__all__ = ['Quadrangle']


class Quadrangle(Polygon):
    """
    Create a patch representing a latitude-longitude quadrangle.
    The edges of the quadrangle lie on two lines of constant longitude and two
    lines of constant latitude (or the equivalent component names in the
    coordinate frame of interest, such as right ascension and declination).
    Note that lines of constant latitude are not great circles.
    Unlike `matplotlib.patches.Rectangle`, the edges of this patch will render
    as curved lines if appropriate for the WCS transformation.
    Parameters
    ----------
    anchor : tuple or `~astropy.units.Quantity`
        This can be either a tuple of two `~astropy.units.Quantity` objects, or
        a single `~astropy.units.Quantity` array with two elements.
    width : `~astropy.units.Quantity`
        The width of the quadrangle in longitude (or, e.g., right ascension)
    height : `~astropy.units.Quantity`
        The height of the quadrangle in latitude (or, e.g., declination)
    resolution : int, optional
        The number of points that make up each side of the quadrangle -
        increase this to get a smoother quadrangle.
    vertex_unit : `~astropy.units.Unit`
        The units in which the resulting polygon should be defined - this
        should match the unit that the transformation (e.g. the WCS
        transformation) expects as input.
    Notes
    -----
    Additional keyword arguments are passed to `~matplotlib.patches.Polygon`
    """

    def __init__(self, anchor, width, height, resolution=100, vertex_unit=u.degree, **kwargs):

        # Extract longitude/latitude, either from a tuple of two quantities, or
        # a single 2-element Quantity.
        longitude, latitude = u.Quantity(anchor).to_value(vertex_unit)

        # Convert the quadrangle dimensions to the appropriate units
        width = width.to_value(vertex_unit)
        height = height.to_value(vertex_unit)

        # Create progressions in longitude and latitude
        lon_seq = longitude + np.linspace(0, width, resolution + 1)
        lat_seq = latitude + np.linspace(0, height, resolution + 1)

        # Trace the path of the quadrangle
        lon = np.concatenate([lon_seq[:-1],
                              np.repeat(lon_seq[-1], resolution),
                              np.flip(lon_seq[1:]),
                              np.repeat(lon_seq[0], resolution)])
        lat = np.concatenate([np.repeat(lat_seq[0], resolution),
                              lat_seq[:-1],
                              np.repeat(lat_seq[-1], resolution),
                              np.flip(lat_seq[1:])])

        # Create polygon vertices
        vertices = np.array([lon, lat]).transpose()

        super().__init__(vertices, **kwargs)
