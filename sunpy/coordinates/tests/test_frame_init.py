# -*- coding: utf-8 -*-

import pytest

import astropy.units as u

from astropy.coordinates import CylindricalRepresentation, UnitSphericalRepresentation, SphericalRepresentation

from .. frames import HelioProjective, HelioGraphicStonyhurst, HelioCentric, HelioGraphicCarrington
from .. representation import UnitSphericalWrap180Representation, SphericalWrap180Representation


def init_frame(frame, args, kwargs):
    if args and kwargs:
        return frame(*args, **kwargs)
    elif args:
        return frame(*args)
    elif kwargs:
        return frame(**kwargs)


"""
These are common 2D params, kwargs are frame specific
"""
two_D_parameters = [([0*u.deg, 0*u.arcsec], None),
                    ([0*u.deg, 0*u.arcsec], {'dateobs':'2011/01/01T00:00:00'}),
                    ([0*u.deg, 0*u.arcsec], {'representation':'unitsphericalwrap180'}),
                    ([0*u.deg, 0*u.arcsec], {'representation':'unitspherical'}),
                    ([UnitSphericalWrap180Representation(0*u.deg, 0*u.arcsec)], None),
                    ([UnitSphericalRepresentation(0*u.deg, 0*u.arcsec)], None),
                    ([UnitSphericalWrap180Representation(0*u.deg, 0*u.arcsec)], {'dateobs':'2011/01/01T00:00:00'})
                   ]
"""
These are common 3D params, kwargs are frame specific
"""
three_D_parameters = [([0*u.deg, 0*u.arcsec, 1*u.Mm], None),
                      ([0*u.deg, 0*u.arcsec, 1*u.Mm], {'dateobs':'2011/01/01T00:00:00'}),
                      ([0*u.deg, 0*u.arcsec, 1*u.Mm], {'representation':'sphericalwrap180'}),
                      ([0*u.deg, 0*u.arcsec, 1*u.Mm], {'representation':'spherical'}),
                      ([SphericalWrap180Representation(0*u.deg, 0*u.arcsec, 1*u.Mm)], None),
                      ([SphericalRepresentation(0*u.deg, 0*u.arcsec, 1*u.Mm)], None),
                      ([SphericalWrap180Representation(0*u.deg, 0*u.arcsec, 1*u.Mm)], {'dateobs':'2011/01/01T00:00:00'})
                     ]

#==============================================================================
### Helioprojective Tests
#==============================================================================

@pytest.mark.parametrize('args, kwargs', two_D_parameters + [(None, {'Tx':0*u.deg, 'Ty':0*u.arcsec})])
def test_create_hpc_2d(args, kwargs):
    hpc1 = init_frame(HelioProjective, args, kwargs)

    # Check we have the right class!
    assert isinstance(hpc1, HelioProjective)
    rep_kwarg = kwargs.get('representation', None) if kwargs else None

    if rep_kwarg and rep_kwarg == 'unitspherical':
        # Check that we have a unitspherical representation
        assert isinstance(hpc1._data, UnitSphericalRepresentation)
    else:
        # Check that we have a 2D wrap180 representation
        assert isinstance(hpc1._data, UnitSphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc1.Tx == 0*u.arcsec
    assert hpc1.Ty == 0*u.arcsec

    # Check the attrs are in the correct default units
    assert hpc1.Tx.unit is u.arcsec
    assert hpc1.Ty.unit is u.arcsec

@pytest.mark.parametrize('args, kwargs', three_D_parameters + [(None, {'Tx':0*u.deg, 'Ty':0*u.arcsec, 'distance':1*u.Mm}),
                                                               ([0*u.deg, 0*u.arcsec], {'distance':1*u.Mm})])
def test_create_3d(args, kwargs):
    hpc1 = init_frame(HelioProjective, args, kwargs)

    # Check we have the right class!
    assert isinstance(hpc1, HelioProjective)
    rep_kwarg = kwargs.get('representation', None) if kwargs else None

    if rep_kwarg and rep_kwarg == 'spherical':
        # Check that we have a unitspherical representation
        assert isinstance(hpc1._data, SphericalRepresentation)
    else:
        # Check that we have a 2D wrap180 representation
        assert isinstance(hpc1._data, SphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc1.Tx == 0*u.arcsec
    assert hpc1.Ty == 0*u.arcsec
    assert hpc1.distance == 1*u.Mm

    # Check the attrs are in the correct default units
    assert hpc1.Tx.unit is u.arcsec
    assert hpc1.Ty.unit is u.arcsec
    assert hpc1.distance.unit is u.km

cylindrical_parameters = [([100*u.km, 25*u.deg, 1*u.Mm], {'representation':'cylindrical'}),
                          ([100*u.km, 25*u.deg, 1*u.Mm], {'dateobs':'2011/01/01T00:00:00',
                                                         'representation':'cylindrical'}),
#                         ([100*u.km, 25*u.deg], {'z':1*u.Mm, 'representation':'cylindrical'}),
#                         (None, {'rho':100*u.km, 'phi':25*u.deg, 'z':1*u.Mm, 'representation':'cylindrical'}),
                         ([CylindricalRepresentation(100*u.km, 25*u.deg, 1*u.Mm)], {'representation':'cylindrical'}),
                         ([CylindricalRepresentation(100*u.km, 25*u.deg, 1*u.Mm)], {'dateobs':'2011/01/01T00:00:00',
                                                                                   'representation':'cylindrical'})
                        ]

@pytest.mark.parametrize('args, kwargs', cylindrical_parameters)
def test_create_cylindrical(args, kwargs):
    hpc1 = init_frame(HelioProjective, args, kwargs)

    # Check we have the right class!
    assert isinstance(hpc1, HelioProjective)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hpc1._data, CylindricalRepresentation)

    # Check the attrs are correct
    assert hpc1.rho == 100*u.km
    assert hpc1.psi == 25*u.deg
    assert hpc1.distance == 1*u.Mm

    # Check the attrs are in the correct default units
    assert hpc1.rho.unit is u.km
    assert hpc1.psi.unit is u.arcsec
    assert hpc1.distance.unit is u.km


#==============================================================================
### Heliographic Stonyhurst Tests
#==============================================================================

@pytest.mark.parametrize('frame', [HelioGraphicStonyhurst, HelioGraphicCarrington])
@pytest.mark.parametrize("args, kwargs", two_D_parameters[:2] + two_D_parameters[4:] + \
                                         [(None, {'lat':0*u.deg, 'lon':0*u.arcsec})])
def test_create_hgs_2d(frame, args, kwargs):
    hgs1 = init_frame(frame, args, kwargs)

    # Check we have the right class!
    assert isinstance(hgs1, frame)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hgs1._data, SphericalWrap180Representation)

    # Check the attrs are correct
    assert hgs1.lon == 0*u.deg
    assert hgs1.lat == 0*u.deg
    assert hgs1.rad == hgs1.RSun

    # Check the attrs are in the correct default units
    assert hgs1.lon.unit is u.deg
    assert hgs1.lat.unit is u.deg
    assert hgs1.rad.unit is u.km


@pytest.mark.parametrize('frame', [HelioGraphicStonyhurst, HelioGraphicCarrington])
@pytest.mark.parametrize("args, kwargs", two_D_parameters[2:4])
def test_create_hgs_force_2d(frame, args, kwargs):
    hgs1 = init_frame(frame, args, kwargs)

    # Check we have the right class!
    assert isinstance(hgs1, frame)

    rep_kwarg = kwargs.get('representation', None) if kwargs else None

    if rep_kwarg == 'unitsphericalwrap180':
        assert isinstance(hgs1._data, UnitSphericalWrap180Representation)
    elif rep_kwarg == 'unitspherical':
        assert isinstance(hgs1._data, UnitSphericalRepresentation)

    # Check the attrs are correct
    assert hgs1.lon == 0*u.deg
    assert hgs1.lat == 0*u.deg

    # Check the attrs are in the correct default units
    assert hgs1.lon.unit is u.deg
    assert hgs1.lat.unit is u.deg

@pytest.mark.parametrize('frame', [HelioGraphicStonyhurst, HelioGraphicCarrington])
@pytest.mark.parametrize("args, kwargs", three_D_parameters+ [(None, {'lat':0*u.deg, 'lon':0*u.arcsec, 'rad':1*u.Mm}),
                                                              ([0*u.deg, 0*u.arcsec], {'rad':1*u.Mm})])
def test_create_hgs_3d(frame, args, kwargs):
    hgs1 = init_frame(frame, args, kwargs)

    # Check we have the right class!
    assert isinstance(hgs1, frame)

    rep_kwarg = kwargs.get('representation', None) if kwargs else None

    if rep_kwarg == 'spherical':
        assert isinstance(hgs1._data, SphericalRepresentation)
    else:
        assert isinstance(hgs1._data, SphericalWrap180Representation)

    # Check the attrs are correct
    assert hgs1.lon == 0*u.deg
    assert hgs1.lat == 0*u.deg
    assert hgs1.rad == 1*u.Mm

    # Check the attrs are in the correct default units
    assert hgs1.lon.unit is u.deg
    assert hgs1.lat.unit is u.deg
    assert hgs1.rad.unit is u.Mm