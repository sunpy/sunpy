"""
This module provides routines for applying solar rotation functions to
mapcubes.
"""

import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u

# SunPy imports
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.image.coalignment import apply_shifts

__author__ = 'J. Ireland'

__all__ = ['calculate_solar_rotate_shift', 'mapcube_solar_derotate']


def calculate_solar_rotate_shift(mc, layer_index=0, **kwargs):
    """
    Calculate the shift that must be applied to each layer
    of a mapcube in order to compensate for solar rotation.  The center
    of the map is used to calculate the position of each mapcube layer.
    Shifts are calculated relative to a specified layer in the mapcube.
    When using this functionality, it is a good idea to check that the
    shifts that were applied to were reasonable and expected.  One way of
    checking this is to animate the original mapcube, animate the derotated
    mapcube, and compare the differences you see to the calculated shifts.

    An example use is as follows.  If you select data from the SDO cutout
    service, it is common to not use the solar tracking implemented by this
    service.  This is because (at time of writing) the solar tracking
    implemented by that service moves the image by single pixels at a time.
    This is not optimal for many use cases, as it introduces artificial jumps
    in the data.  So with solar tracking not chosen, the selected area is
    like a window through which you can see the Sun rotating underneath.

    Parameters
    ----------
    mc : `sunpy.map.MapCube`
        The input mapcube.
    layer_index : int
        The index layer.  Shifts are calculated relative to the time of
        this layer.
    ``**kwargs``
        These keywords are passed to the function
        `sunpy.physics.differential_rotation.solar_rotate_coordinate`.

    Returns
    -------
    x, y : `~astropy.units.Quantity`, ~astropy.units.Quantity`
        The shifts relative to the index layer that can be applied
        to the input mapcube in order to compensate for solar rotation.
        The shifts are given in arcseconds as understood in helioprojective
        coordinates systems.

    """
    # Size of the data
    nt = len(mc.maps)

    # Storage for the shifts in arcseconds
    xshift_arcseconds = np.zeros(nt) * u.arcsec
    yshift_arcseconds = np.zeros_like(xshift_arcseconds)

    # Layer that
    rotate_to_this_layer = mc.maps[layer_index]

    # Calculate the rotations and the shifts
    for i, m in enumerate(mc):
        # Calculate the rotation of the center of the map 'm' at its
        # observation time to the observation time of the reference layer
        # indicated by "layer_index".
        new_coordinate = solar_rotate_coordinate(m.center,
                                                 rotate_to_this_layer.date,
                                                 new_observer_location=rotate_to_this_layer.observer_coordinate,
                                                 **kwargs)

        # Calculate the shift in arcseconds
        xshift_arcseconds[i] = new_coordinate.Tx - rotate_to_this_layer.center.Tx
        yshift_arcseconds[i] = new_coordinate.Ty - rotate_to_this_layer.center.Ty

    return {"x": xshift_arcseconds, "y": yshift_arcseconds}


def mapcube_solar_derotate(mc, layer_index=0, clip=True, shift=None, **kwargs):
    """
    Move the layers in a mapcube according to the input shifts.
    If an input shift is not given, the shifts due to
    solar rotation relative to an index layer is calculated and
    applied.  When using this functionality, it is a good idea to check
    that the shifts that were applied to were reasonable and expected.
    One way of checking this is to animate the original mapcube, animate
    the derotated mapcube, and compare the differences you see to the
    calculated shifts.

    Parameters
    ----------
    mc : `sunpy.map.MapCube`
        A mapcube of shape (ny, nx, nt), where nt is the number of layers in
        the mapcube.

    layer_index : int
        Solar derotation shifts of all maps in the mapcube are assumed
        to be relative to the layer in the mapcube indexed by layer_index.

    clip : bool
        If True, then clip off x, y edges in the datacube that are potentially
        affected by edges effects.

    ``**kwargs``
        These keywords are passed to the function
        `sunpy.physics.solar_rotation.calculate_solar_rotate_shift`.

    Returns
    -------
    output : `sunpy.map.MapCube`
        The results of the shifts applied to the input mapcube.

    Examples
    --------
    >>> from sunpy.physics.solar_rotation import mapcube_solar_derotate
    >>> import sunpy.data.sample
    >>> map1 = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> map2 = sunpy.map.Map(sunpy.data.sample.EIT_195_IMAGE)
    >>> mc = sunpy.map.Map([map1, map2], cube=True)
    >>> derotated_mc = mapcube_solar_derotate(mc)
    >>> derotated_mc = mapcube_solar_derotate(mc, layer_index=-1)
    >>> derotated_mc = mapcube_solar_derotate(mc, clip=False)
    """

    # Size of the data
    nt = len(mc.maps)

    # Storage for the pixel shifts and the shifts in arcseconds
    xshift_keep = np.zeros((nt)) * u.pix
    yshift_keep = np.zeros_like(xshift_keep)

    # If no shifts are passed in, calculate them.  Otherwise,
    # use the shifts passed in.
    if shift is None:
        shift = calculate_solar_rotate_shift(mc, layer_index=layer_index, **kwargs)
    xshift_arcseconds = shift['x']
    yshift_arcseconds = shift['y']

    # Calculate the pixel shifts
    for i, m in enumerate(mc):
        xshift_keep[i] = xshift_arcseconds[i] / m.scale[0]
        yshift_keep[i] = yshift_arcseconds[i] / m.scale[1]

    # Apply the pixel shifts and return the mapcube
    return apply_shifts(mc, yshift_keep, xshift_keep, clip=clip)
