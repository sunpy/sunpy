"""
This module provides routines for applying solar differential functions to
mapcubes.
"""

import numpy as np
from astropy import units as u

# SunPy imports
from sunpy.physics.transforms.differential_rotation import rot_hpc
from sunpy.image.coalignment import apply_shifts

__author__ = 'J. Ireland'

__all__ = ['mapcube_solar_derotate']


def mapcube_solar_derotate(mc, layer_index=0, clip=True,
                               return_displacements_only=False,
                               with_displacements=False, **kwargs):
    """Move the layers in a mapcube according to solar rotation.  The center
    of the map is used to calculate the position of each mapcube layer.  Shifts
    are calculated relative to a specified layer in the mapcube.
    When using this functionality, it is a good idea to check that the
    shifts that were applied to were reasonable and expected.  One way of
    checking this is to animate the original mapcube, animate the derotated
    mapcube, and compare the differences you see to the calculated shifts.


    Parameters
    ----------
    mc : sunpy.map.MapCube
        A mapcube of shape (ny, nx, nt), where nt is the number of layers in
        the mapcube.

    layer_index : int
        Solar derotation displacements of all maps in the mapcube are assumed
        to be relative to the layer in the mapcube indexed by layer_index.

    clip : bool
        If True, then clip off x, y edges in the datacube that are potentially
        affected by edges effects.

    return_displacements_only : bool
        If True return ONLY the x and y displacements applied to the input
        data in units of arcseconds.  The return value is a dictionary of the
        form {"x": xdisplacement, "y": ydisplacement}.

    with_displacements : bool
        If True, return the x and y displacements applied to the input data in
        the same format as that returned using the return_displacements_only
        option, along with the derotated mapcube.  The format of the return is
        (mapcube, displacements).

    Returns
    -------
    output : {sunpy.map.MapCube | dict | tuple}
        The results of the mapcube coalignment.  The output depends on the
        value of the parameters "return_displacements_only" and
        "with_displacements".
"""
    # Size of the data
    ny = mc.maps[layer_index].shape[0]
    nx = mc.maps[layer_index].shape[1]
    nt = len(mc.maps)

    # Storage for the pixel shifts and the shifts in arcseconds
    xshift_keep = np.zeros((nt))
    yshift_keep = np.zeros_like(xshift_keep)

    # Calculate the rotations and the shifts
    for i, m in enumerate(mc):
        newx, newy = rot_hpc(xc, yc, mc.maps[layer_index].date, m.date, **kwargs)

    # Return only the displacements
    if return_displacements_only:
        return {"x": xshift_arcseconds * u.arcsec, "y": yshift_arcseconds * u.arcsec}

    # Apply the pixel shifts
    newmc = apply_shifts(mc, -yshift_keep * u.pix, -xshift_keep * u.pix, clip=clip)

   # Return the mapcube, or optionally, the mapcube and the displacements
    # used to create the mapcube.
    if with_displacements:
        return newmc, {"x": xshift_arcseconds * u.arcsec, "y": yshift_arcseconds * u.arcsec}
    else:
        return newmc
