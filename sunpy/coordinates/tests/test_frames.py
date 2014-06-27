"""
Frame class tests for the coordinates framework.
Modelled after tests in Astropy.
@author: Pritish C. (VaticanCameos)
"""

# NumPy
import numpy as np
from numpy import testing as npt

# Astropy
from astropy import units as u
from astropy.coordinates import (SphericalRepresentation,
                                 CylindricalRepresentation,
                                 CartesianRepresentation,
                                 UnitSphericalRepresentation)

# Pytest
import pytest

# SunPy
from ..frames import (HelioGraphicStonyhurst, HelioGraphicCarrington,
                      HelioCentric, HelioProjective)

def test_frame_attributes():
    from ..frames import HelioProjective
    from astropy.coordinates.baseframe import FrameAttribute

    class MyProj(HelioProjective):
        # Inherit D0, change d to something else, create a newattr.
        d = FrameAttribute(default=(2*u.au).to(u.km))
        newattr = FrameAttribute(default='anattr')

    myproj = MyProj()
    assert myproj.D0 == (1*u.au).to(u.km)
    assert myproj.d == (2*u.au).to(u.km)
    assert myproj.newattr == 'anattr'
    assert set(myproj.get_frame_attr_names()) == set(['D0',
                                                      'd',
                                                      'newattr'])

    myproj = MyProj(D0=1*u.km, d=2*u.km, newattr='changed')
    assert myproj.D0 == 1*u.km
    assert myproj.d == 2*u.km
    assert myproj.newattr == 'changed'

    
def test_create_spherical_data_frames():
    for frame in [HelioGraphicStonyhurst, HelioGraphicCarrington,
                  HelioProjective]:

        # Using representations
        coord1 = frame(SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc))

        # Using preferred names
        coord2 = frame(1*u.deg, 2*u.deg, 3*u.kpc)

        assert coord1.data.lat == coord2.data.lat
        assert coord1.data.lon == coord2.data.lon
        assert coord1.data.distance == coord2.data.distance

        # Check as properties
        if isinstance(coord1, HelioGraphicStonyhurst):
            npt.assert_allclose(coord1.hlon, coord2.hlon)
            npt.assert_allclose(coord1.hlat, coord2.hlat)
            npt.assert_allclose(coord1.rad, coord2.rad)

            with pytest.raises(AttributeError):
                coord1.hlat = [10.]*u.deg
        else:
            npt.assert_allclose(coord1.Tx, coord2.Tx)
            npt.assert_allclose(coord1.Ty, coord2.Ty)
            npt.assert_allclose(coord1.zeta, coord2.zeta)

            with pytest.raises(AttributeError):
                coord1.Tx = 5*u.deg

        
def test_create_data_cartesian_frames():
    # Using representations
    coord1 = HelioCentric(CartesianRepresentation(1*u.km, 2*u.km, 3*u.km))

    # Using frames
    coord2 = HelioCentric(x=1*u.km, y=2*u.km, z=3*u.km)

    assert coord1.data.x == coord2.data.x
    assert coord1.data.y == coord2.data.y
    assert coord1.data.z == coord2.data.z

    # Check as properties
    npt.assert_allclose(coord1.x, coord2.x)
    npt.assert_allclose(coord1.y, coord2.y)
    npt.assert_allclose(coord1.z, coord2.z)

    with pytest.raises(AttributeError):
        coord1.x = 1*u.km
