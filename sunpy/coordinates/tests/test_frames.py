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
from astropy.coordinates import SkyCoord

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
                                                      'newattr',
                                                      'dateobs',
                                                      'L0','B0'])

    myproj = MyProj(D0=1*u.km, d=2*u.km, newattr='changed')
    assert myproj.D0 == 1*u.km
    assert myproj.d == 2*u.km
    assert myproj.L0 == 0*u.deg
    assert myproj.B0 == 0*u.deg
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
    assert len(hgs.get_frame_attr_names()) == 3
    # Heliographic frames are either completely empty (with the dateobs kwarg),
    # or true 3D.
    # Nothing in between.

    hcc = HelioCentric()
    assert hcc.D0 == HelioCentric.get_frame_attr_names()['D0']
    assert len(hcc.get_frame_attr_names()) == 4

    hp = HelioProjective()
    assert hp.D0 == HelioProjective.get_frame_attr_names()['D0']
    assert len(hp.get_frame_attr_names()) == 4

def test_frame_repr():
    # Tests the repr() of a frame.

    hgc = HelioGraphicCarrington()
    assert repr(hgc).startswith('<HelioGraphicCarrington Frame: ')

    hcc = HelioCentric()
    assert repr(hcc).startswith('<HelioCentric Frame: ')

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

    hgs = SkyCoord(HelioGraphicStonyhurst(hlon=[1,2]*u.deg, hlat=[3,4]*u.deg,
                                 dateobs="2011/01/01T00:00:45"))
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

@pytest.mark.parametrize("hgs1, hgs2", [
    (HelioGraphicStonyhurst(1*u.deg, 1*u.deg, dateobs=datetime.now()),
     HelioGraphicStonyhurst(1*u.deg, 1*u.deg, 1*u.km, dateobs=datetime.now())),
    (HelioGraphicStonyhurst(SphericalWrap180Representation(1*u.deg, 1*u.deg, 1*u.km), dateobs=datetime.now()),
     HelioGraphicStonyhurst(SphericalWrap180Representation(1*u.deg, 1*u.deg, 1*u.km))),
    (HelioGraphicStonyhurst(1*u.deg, hlat=1*u.deg, rad=1*u.km, dateobs=datetime.now()),
     HelioGraphicStonyhurst(1*u.deg, hlat=1*u.deg, rad=1*u.km)),
    (HelioGraphicStonyhurst(hlon=1*u.deg, hlat=1*u.deg, rad=1*u.km, dateobs=datetime.now()),
     HelioGraphicStonyhurst(hlon=1*u.deg, hlat=1*u.deg, rad=1*u.km)),
    (HelioGraphicStonyhurst(dateobs=datetime.now()), None)])
def test_heliographic_stonyhurst(hgs1, hgs2):
    # A test function for testing the several ways
    # to create HGS frames.

    if hgs2 is not None:
        npt.assert_allclose(hgs1.hlon, hgs2.hlon)
        npt.assert_allclose(hgs1.hlat, hgs2.hlat)
        assert hgs1.dateobs != hgs2.dateobs
    else:
        # Deals with empty frame case.
        assert len(hgs1.get_frame_attr_names()) == 3 # Should only be dateobs
        assert hgs1.dateobs is not None
        

def test_transform_accuracy():
    """
    This tests the accuracy of transforms.
    Some input values are fed and output is checked to be
    equal within an allowed tolerance level to expected output.
    Here expected output is the output that was obtained by
    working out the equations by hand.
    """
    from sunpy import sun as s
    from sunpy.coordinates.representation import Longitude180
    from numpy import sqrt, arcsin, arctan
    import numpy
    
    RSun = s.constants.constant('radius').si.to(u.km)
    DSun = s.constants.constant('mean distance').si.to(u.km)
    diff = (1*u.au).to(u.km) - RSun
    
    # First work on (0,0,RSun) case.    
    # Explicitly pass RSun because we're not testing the constructor's
    # defaulting capability here.
    sc_zero = SkyCoord(0*u.deg, 0*u.deg, RSun, frame='heliographicstonyhurst',
                       dateobs='2011/01/01T00:00:45')
    dateobs_calc = s.heliographic_solar_center(sc_zero.dateobs)
    sc_zero_hp = sc_zero.transform_to('helioprojective')
    sc_zero_hcc = sc_zero.transform_to('heliocentric')
    sc_zero_hgc = sc_zero.transform_to('heliographiccarrington')
    
    npt.assert_allclose(sc_zero.hlon.value, sc_zero_hp.Tx.value)
    npt.assert_allclose(sc_zero.hlat.value, sc_zero_hp.Ty.value)
    npt.assert_allclose(sc_zero_hp.distance, diff)  

    npt.assert_allclose(sc_zero.hlon.value, sc_zero_hcc.x.value)
    npt.assert_allclose(sc_zero.hlat.value, sc_zero_hcc.y.value)
    npt.assert_allclose(sc_zero_hcc.z, RSun)

    npt.assert_allclose(sc_zero.hlon, sc_zero_hgc.hlon - Longitude180(dateobs_calc[0]))
    npt.assert_allclose(sc_zero.hlat, sc_zero_hgc.hlat)
    npt.assert_allclose(sc_zero.rad, sc_zero_hgc.rad)
    
    # Cadair's analytical test
    # See: http://nbviewer.ipython.org/gist/Cadair/63d405b956c3478bfa64
    # Also, see: http://nbviewer.ipython.org/gist/Cadair/d795984b877f16c8eeb2
    # Params: lon=45deg, lat=45deg, rad=RSun, L0=0deg, B0=0deg
    
    sc_hgs = SkyCoord(45*u.deg, 45*u.deg, RSun, L0=0*u.deg, B0=0*u.deg,
                      frame='heliographicstonyhurst',
                      dateobs='2011/01/01T00:00:45')
                      
    # Now in HCC, we should have x=RSun/2, y=RSun/root2, z=RSun/2.
    
    expect_x = RSun/2
    expect_y = RSun/sqrt(2)
    expect_z = RSun/2
    
    sc_hcc = sc_hgs.transform_to('heliocentric')
    
    npt.assert_allclose(sc_hcc.x, expect_x)
    npt.assert_allclose(sc_hcc.y, expect_y)
    npt.assert_allclose(sc_hcc.z, expect_z)
    
    # Now in HPC, we have a complex bunch of equations. 
    
    expect_Ty = arcsin((RSun/sqrt(2))/(np.sqrt(((3*RSun**2)/4) + (DSun - (RSun/2))**2))).to(u.deg)
    expect_Tx = arctan(1/(2*(DSun/RSun)-1)).to(u.deg)
    expect_d = np.sqrt(((3*RSun**2)/4) + (DSun - (RSun/2))**2).to(u.km)
    
    sc_hpc1 = sc_hgs.transform_to('helioprojective')
    sc_hpc2 = sc_hcc.transform_to('helioprojective')
    
    for coord in [sc_hpc1, sc_hpc2]:
        npt.assert_allclose(coord.Tx.to(u.deg), expect_Tx)
        npt.assert_allclose(coord.distance, expect_d)
        npt.assert_allclose(coord.Ty.to(u.deg), expect_Ty)
        
from sunpy import sun as s
RSun = s.constants.radius.si.to(u.km)   

@pytest.mark.parametrize("input, expected, extras, to",
                         [([40.0, 32.0] * u.arcsec, [28748691, 22998953] * u.m, 
                           {'frame': 'helioprojective'}, 'heliocentric'),
                          ([40.0, 32.0] * u.arcsec, [28748691, 22998953] * u.m,
                           {'distance': 0.5 * u.au, 'frame': 'helioprojective'},
                            'heliocentric'),
                          ([28748691, 22998953, 0] * u.m, [40.0, 32.0] * u.arcsec,
                           {'frame': 'heliocentric'}, 'helioprojective'),
                          ([13.0, 58.0, 0] * u.m, [1.0791282e-06*u.deg, -7.0640732*u.deg,
                           RSun],
                           {'frame': 'heliocentric'}, 'heliographicstonyhurst')])
def test_wcs_numbers(input, expected, extras, to):
    dateobs = '2011/01/01T00:00:45'
    extras['dateobs'] = dateobs
    rtol = 1e-10
    sc = SkyCoord(*input, **extras)
    
    sc_trans = sc.transform_to(to)
    
    if sc_trans.representation is SphericalWrap180Representation:
        npt.assert_allclose(sc_trans.spherical.lon, expected[0].to(u.deg))
        npt.assert_allclose(sc_trans.spherical.lat, expected[1].to(u.deg))
        if expected[2] is not None:
            npt.assert_allclose(sc_trans.spherical.distance, expected[2].to(u.km))
    elif sc_trans.representation is CartesianRepresentation:
        npt.assert_allclose(sc_trans.cartesian.x, expected[0].to(u.km))
        npt.assert_allclose(sc_trans.cartesian.y, expected[1].to(u.km))
        if expected[2] is not None:
            npt.assert_allclose(sc_trans.cartesian.z, expected[2].to(u.km))