# -*- coding: utf-8 -*-

import astropy.units as u

from astropy.coordinates import CylindricalRepresentation

from .. frames import HelioProjective, HelioGraphicStonyhurst, HelioCentric
from .. representation import UnitSphericalWrap180Representation, SphericalWrap180Representation

def test_create_hpc_2d():
    hpc1 = HelioProjective(0.1*u.deg, 200*u.arcsec)

    # Check we have the right class!
    assert isinstance(hpc1, HelioProjective)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hpc1._data, UnitSphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc1.Tx == 360*u.arcsec
    assert hpc1.Ty == 200*u.arcsec

    # Check the attrs are in the correct default units
    assert hpc1.Tx.unit is u.arcsec
    assert hpc1.Ty.unit is u.arcsec

def test_create_3d():
    hpc1 = HelioProjective(100*u.arcsec, 200*u.arcsec, 1*u.Mm)

    # Check we have the right class!
    assert isinstance(hpc1, HelioProjective)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hpc1._data, SphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc1.Tx == 100*u.arcsec
    assert hpc1.Ty == 200*u.arcsec
    assert hpc1.distance == 1*u.Mm

    # Check the attrs are in the correct default units
    assert hpc1.Tx.unit is u.arcsec
    assert hpc1.Ty.unit is u.arcsec
    assert hpc1.distance.unit is u.km


def test_create_cylindrical():
    hpc1 = HelioProjective(100*u.km, 25*u.deg, 1*u.Mm,
                           representation='cylindrical')

    # Check we have the right class!
    assert isinstance(hpc1, HelioProjective)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hpc1._data, CylindricalRepresentation)

    # Check the attrs are correct
    assert hpc1.Trho == 100*u.km
    assert hpc1.psi == 25*u.deg
    assert hpc1.distance == 1*u.Mm

    # Check the attrs are in the correct default units
    assert hpc1.Trho.unit is u.km
    assert hpc1.psi.unit is u.arcsec
    assert hpc1.distance.unit is u.km
