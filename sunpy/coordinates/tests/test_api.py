"""
Test cases for the upcoming cooordinates API.
Tests will pass gradually as features are
added.
@Author: Pritish C.(VaticanCameos)
"""

import numpy as np
from numpy import testing as npt
from astropy.extern import six
import pytest

from astropy import units as u

def test_frame_api():
    from astropy.coordinates import SphericalRepresentation, \
         CylindricalRepresentation, CartesianRepresentation
    from ..sunpy_frames import HelioGraphic, HelioCentric, \
         HelioProjective
    # This method will test the use of the Coordinates API
    # in a way similar to that of Astropy's testing methods.

    hgic = HelioGraphic(SphericalRepresentation(lon=3*u.deg, lat=10*u.deg))

    # Check the data attribute.
    assert hgic.data.lon == 3*u.deg
    assert hgic.data.lat == 10*u.deg

    # Representation objects give us the actual positional information.
    npt.assert_allclose(hgic.represent_as(SphericalRepresentation).lat.to(u.deg), 10*u.deg)
    npt.assert_allclose(hgic.spherical.lat.to(u.deg), 10*u.deg)
    assert hgic.represent_as(CartesianRepresentation).z.value > 0

    # Now assert the same with hgic's attributes and not the representation's.
    npt.assert_allclose(hgic.lat.to(u.deg), 10*u.deg)
    npt.assert_allclose(hgic.lon.to(u.deg), 3*u.deg)

    assert hgic.representation is SphericalRepresentation

    # Now we try a low-level initialization of hgic.
    hgic_2 = HelioGraphic(lon=10*u.deg, lat=8*u.deg, radius=1*u.kpc)
    npt.assert_allclose(hgic.lon.to(u.deg), hgic_2.lon.to(u.deg))

    # Try calculating the separation and testing against it.
    # This is the on-sky separation.
    coord1 = HelioGraphic(lon=0*u.deg, lat=0*u.deg)
    coord2 = HelioGraphic(lon=0*u.deg, lat=1*u.deg)
    assert coord1.separation(coord2).degree == 1.0

    # Now try the same with the 3D separation.
    coord3 = HelioGraphic(lon=0*u.deg, lat=0*u.deg, radius=1*u.kpc)
    coord4 = HelioGraphic(lon=0*u.deg, lat=0*u.deg, radius=2*u.kpc)
    assert coord3.separation_3d(coord4).kpc == 1.0

