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
from astropy.coordinates.representation import (SphericalRepresentation,
                                                CylindricalRepresentation,
                                                CartesianRepresentation)
from ..representation import SphericalWrap180Representation

def test_frame_api():
    from ..frames import (HelioGraphicStonyhurst, HelioCentric,
                          HelioProjective)
    # This method will test the use of the Coordinates API
    # in a way similar to that of Astropy's testing methods.

    hgic = HelioGraphicStonyhurst(SphericalRepresentation(lon=3*u.deg, lat=10*u.deg,
                                                          distance=1*u.km))

    # Check the data attribute.
    assert hgic.data.lon == 3*u.deg
    assert hgic.data.lat == 10*u.deg

    # Representation objects give us the actual positional information.
    npt.assert_allclose(hgic.represent_as(SphericalRepresentation).lat.to(u.deg), 10*u.deg)
    npt.assert_allclose(hgic.spherical.lat.to(u.deg), 10*u.deg)
    assert hgic.represent_as(CartesianRepresentation).z.value > 0

    # Now assert the same with hgic's attributes and not the representation's.
    npt.assert_allclose(hgic.hlat.to(u.deg), 10*u.deg)
    npt.assert_allclose(hgic.hlon.to(u.deg), 3*u.deg)

    assert hgic.representation is SphericalWrap180Representation

    # Now we try a low-level initialization of hgic.
    hgic_2 = HelioGraphicStonyhurst(3*u.deg, 8*u.deg, 1*u.kpc)
    npt.assert_allclose(hgic.hlon.to(u.deg), hgic_2.hlon.to(u.deg))

    # Try calculating the separation and testing against it.
    # This is the on-sky separation.
    # The radius attribute, if not specified, is defaulted
    # to the solar radius.
    #coord1 = HelioGraphicStonyhurst(hlon=0*u.deg, hlat=0*u.deg)
    #coord2 = HelioGraphicStonyhurst(hlon=0*u.deg, hlat=1*u.deg)
    #assert coord1.separation(coord2).degree == 1.0

    # Now try the same with the 3D separation.
    coord3 = HelioGraphicStonyhurst(rad=1*u.kpc, hlon=1*u.deg, hlat=1*u.deg)
    coord4 = HelioGraphicStonyhurst(rad=2*u.kpc, hlon=1*u.deg, hlat=1*u.deg)
    assert coord3.separation_3d(coord4).kpc == 1.0

def test_highlevel_api():
    from astropy.coordinates import SkyCoord
    from ..frames import HelioGraphicStonyhurst, HelioCentric
    
    # This method tests the high-level API as is done in Astropy.
    # SkyCoord should be SunPy-ready by the time the coordinates
    # API is in the codebase.

    sc = SkyCoord(SphericalRepresentation(lon=10*u.deg, lat=10*u.deg,
                                          distance=1*u.kpc),
                  frame="heliographicstonyhurst",dateobs="2011/01/01T00:00:45")

    # Heliocentric coordinates are in kilometres.
    sc = SkyCoord(hlon=10*u.deg, hlat=10*u.deg, frame="heliographicstonyhurst")
    sc = SkyCoord(x=10*u.km, y=10*u.km, z=10*u.km, frame="heliocentric")

    # One can initialize using low-level objects.
    sc = SkyCoord(HelioGraphicStonyhurst(hlon=8*u.deg, hlat=10*u.deg, dateobs=
    "2011/01/01T00:00:45"))

    # An error is induced as a high-level object needs position data.
    # Frames can be initialized without this data, SkyCoord cannot.
    with pytest.raises(ValueError):
        sc = SkyCoord(frame="heliocentric")

    # The underlying frame object of the high-level object, when
    # accessed in a call to `repr`, is printed in the following
    # way - '<HelioGraphicStonyhurst Coordinate: lon=10*u.deg, lat=10*u.deg>'
    string = repr(sc.frame)
    assert '<HelioGraphicStonyhurst Coordinate: B0=' in string
    assert 'hlon=' in string
    assert 'deg, hlat=' in string
    assert 'rad=' in string
    assert 'L0=' in string
    assert 'km>' in string

    # Similarly, `repr(sc)` should look like -:
    # '<SkyCoord (HelioGraphicStonyhurst): lon=10*u.deg, lat=10*u.deg>'
    string = repr(sc)
    assert '<SkyCoord (HelioGraphicStonyhurst): B0=' in string
    assert '<HelioGraphicStonyhurst Coordinate: B0=' in string
    assert 'hlon=' in string
    assert 'deg, hlat=' in string
    assert 'rad=' in string
    assert 'L0=' in string
    assert 'km>' in string

    # Transformation between frames is delegated by SkyCoord
    # to lower-level frame classes.

    sc_heliocentric = sc.transform_to("heliocentric")
    assert sc_heliocentric.representation is CartesianRepresentation

    # SkyCoord also contains transformations to other frames as part of
    # its attributes.
    sc = SkyCoord(hlon=8*u.deg, hlat=10*u.deg, frame="heliographicstonyhurst",
                  dateobs="2011/01/01T00:00:45")
    sc_helioprojective = sc.helioprojective
    assert repr(sc_helioprojective).startswith('<SkyCoord (HelioProjective):')
