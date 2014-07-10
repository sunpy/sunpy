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
from astropy.coordinates.representation import (SphericalRepresentation,
                                                CylindricalRepresentation,
                                                CartesianRepresentation,
                                                UnitSphericalRepresentation)

# Pytest
import pytest

# SunPy
from ..frames import (HelioGraphicStonyhurst, HelioGraphicCarrington,
                      HelioCentric, HelioProjective)
from ..representation import SphericalWrap180Representation
from datetime import datetime

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

def test_ordered_data_frames():
    # Need a similar test function for transformations
    # involving heliographic frames - due to L0/B0.

    # Tolerance level.
    TOL = 1e-10*u.deg

    hgs = HelioGraphicStonyhurst(1*u.deg, 2*u.deg)
    assert (hgs.hlon - 1*u.deg) < TOL
    assert (hgs.hlat - 2*u.deg) < TOL

    hgc = HelioGraphicCarrington(1*u.deg, 2*u.deg)
    assert (hgc.hlon - 1*u.deg) < TOL
    assert (hgc.hlat - 2*u.deg) < TOL

    hp = HelioProjective(1*u.deg, 2*u.deg, 3*u.km)
    assert (hp.Tx - 1*u.deg) < TOL
    assert (hp.Ty - 2*u.deg) < TOL

    with pytest.raises(TypeError):
        HelioCentric(1*u.km, 2*u.km, 3*u.km, 4*u.km)

    with pytest.raises(TypeError):
        cr = CartesianRepresentation(1*u.km, 2*u.km, 3*u.km)
        HelioCentric(cr, 4*u.km)

def test_nodata_frames():
    # Tests frames which have no data.

    hgs = HelioGraphicStonyhurst()
    assert len(hgs.get_frame_attr_names()) == 1
    # Heliographic frames are either completely empty (with the dateobs kwarg),
    # or true 3D.
    # Nothing in between.

    hcc = HelioCentric()
    assert hcc.D0 == HelioCentric.get_frame_attr_names()['D0']
    assert len(hcc.get_frame_attr_names()) == 1

    hp = HelioProjective()
    assert hp.D0 == HelioProjective.get_frame_attr_names()['D0']
    assert len(hp.get_frame_attr_names()) == 1

def test_frame_repr():
    # Tests the repr() of a frame.

    hgc = HelioGraphicCarrington()
    assert repr(hgc).startswith('<HelioGraphicCarrington Frame: dateobs=')

    hcc = HelioCentric()
    assert repr(hcc).startswith('<HelioCentric Frame: D0=')

    hp = HelioProjective()
    assert '<HelioProjective Frame:' in repr(hp)
    assert 'D0=' in repr(hp)
    
    hgs_1 = HelioGraphicStonyhurst(1*u.deg, 2*u.deg)
    hgs_2 = HelioGraphicStonyhurst(1*u.deg, 2*u.deg, 3*u.km)

    string1 = repr(hgs_1)
    string2 = repr(hgs_2)
    assert '<HelioGraphicStonyhurst Coordinate: ' in string1
    assert 'dateobs=' in string1
    assert 'hlon=1.0 deg, hlat=2.0 deg, rad=695508.0 km>' in string1
    assert '<HelioGraphicStonyhurst Coordinate: ' in string2
    assert 'dateobs=' in string2
    assert 'hlon=1.0 deg, hlat=2.0 deg, rad=3.0 km>' in string2

def test_realize_frames():
    # Tests for the realize_frame() method.

    rep = SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.km)

    hgs_1 = HelioGraphicStonyhurst()
    hgs_2 = hgs_1.realize_frame(rep)

    assert not hgs_1.has_data
    assert hgs_2.has_data

    hp_1 = HelioProjective(D0=2*u.km)
    hp_2 = hp_1.realize_frame(rep)

    assert not hp_1.has_data
    assert hp_2.has_data

    assert hp_1.D0 == hp_2.D0
    assert hp_2.D0 != HelioProjective.get_frame_attr_names()['D0']

def test_transform_architecture():
    # This does not test the accuracy of transforms.

    hgs = HelioGraphicStonyhurst(hlon=[1,2]*u.deg,
                                 hlat=[3,4]*u.deg)
    hcc = hgs.transform_to(HelioCentric)
    hgs_2 = hcc.transform_to(HelioGraphicStonyhurst)
    hgc = hgs_2.transform_to(HelioGraphicCarrington)

    assert hgs_2.data.__class__ != UnitSphericalRepresentation
    assert hgc.data.__class__ != UnitSphericalRepresentation

    npt.assert_allclose(hgs.hlon, hgs_2.hlon)
    npt.assert_allclose(hgs.hlat, hgs_2.hlat)

    # Self transforms

    hgs_3 = hgs.transform_to(HelioGraphicStonyhurst)

    npt.assert_allclose(hgs.hlon, hgs_3.hlon)
    npt.assert_allclose(hgs.hlat, hgs_3.hlat)
    npt.assert_allclose(hgs.rad, hgs_3.rad)

def test_heliographic_stonyhurst():
    # A test function for testing the several ways
    # to create HGS frames.

    HelioGraphicStonyhurst(1*u.deg, 1*u.deg, dateobs=datetime.now())
    HelioGraphicStonyhurst(1*u.deg, 1*u.deg, 1*u.km, dateobs=datetime.now())
    HelioGraphicStonyhurst(SphericalWrap180Representation(1*u.deg, 1*u.deg, 1*u.km), dateobs=datetime.now())
    HelioGraphicStonyhurst(SphericalWrap180Representation(1*u.deg, 1*u.deg, 1*u.km))
    HelioGraphicStonyhurst(dateobs=datetime.now())
    HelioGraphicStonyhurst(1*u.deg, hlat=1*u.deg, rad=1*u.km, dateobs=datetime.now())
    HelioGraphicStonyhurst(1*u.deg, hlat=1*u.deg, rad=1*u.km)
    HelioGraphicStonyhurst(hlon=1*u.deg, hlat=1*u.deg, rad=1*u.km, dateobs=datetime.now())
    HelioGraphicStonyhurst(hlon=1*u.deg, hlat=1*u.deg, rad=1*u.km)

    

    
