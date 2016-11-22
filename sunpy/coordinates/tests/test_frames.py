# -*- coding: utf-8 -*-

import numpy as np

import pytest

import astropy.units as u

from astropy.tests.helper import assert_quantity_allclose

from astropy.coordinates import (
    CylindricalRepresentation, UnitSphericalRepresentation,
    SphericalRepresentation, CartesianRepresentation)

from ... import sun
from ..frames import Helioprojective, HeliographicStonyhurst, Heliocentric, HeliographicCarrington
from ..representation import UnitSphericalWrap180Representation, SphericalWrap180Representation

RSUN_METERS = sun.constants.get('radius').si.to(u.m)
DSUN_METERS = sun.constants.get('mean distance').si.to(u.m)


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
two_D_parameters = [
    ([0 * u.deg, 0 * u.arcsec], None),
    ([0 * u.deg, 0 * u.arcsec], {'dateobs': '2011/01/01T00:00:00'}),
    ([0 * u.deg, 0 * u.arcsec], {'representation': 'unitsphericalwrap180'}),
    ([0 * u.deg, 0 * u.arcsec], {'representation': 'unitspherical'}),
    ([UnitSphericalWrap180Representation(0 * u.deg, 0 * u.arcsec)], None),
    ([UnitSphericalRepresentation(0 * u.deg, 0 * u.arcsec)], None), (
        [UnitSphericalWrap180Representation(0 * u.deg, 0 * u.arcsec)],
        {'dateobs': '2011/01/01T00:00:00'})
]
"""
These are common 3D params, kwargs are frame specific
"""
three_D_parameters = [
    ([0 * u.deg, 0 * u.arcsec, 1 * u.Mm], None),
    ([0 * u.deg, 0 * u.arcsec, 1 * u.Mm], {'dateobs': '2011/01/01T00:00:00'}),
    ([0 * u.deg, 0 * u.arcsec, 1 * u.Mm], {'representation': 'sphericalwrap180'
                                           }),
    ([0 * u.deg, 0 * u.arcsec, 1 * u.Mm], {'representation': 'spherical'}),
    ([SphericalWrap180Representation(0 * u.deg, 0 * u.arcsec, 1 * u.Mm)],
     None),
    ([SphericalRepresentation(0 * u.deg, 0 * u.arcsec, 1 * u.Mm)], None), (
        [SphericalWrap180Representation(0 * u.deg, 0 * u.arcsec, 1 * u.Mm)],
        {'dateobs': '2011/01/01T00:00:00'})
]

#==============================================================================
# Helioprojective Tests
#==============================================================================


@pytest.mark.parametrize('args, kwargs',
                         two_D_parameters + [(None, {'Tx': 0 * u.deg,
                                                     'Ty': 0 * u.arcsec})])
def test_create_hpc_2d(args, kwargs):
    hpc1 = init_frame(Helioprojective, args, kwargs)

    # Check we have the right class!
    assert isinstance(hpc1, Helioprojective)
    rep_kwarg = kwargs.get('representation', None) if kwargs else None

    if rep_kwarg and rep_kwarg == 'unitspherical':
        # Check that we have a unitspherical representation
        assert isinstance(hpc1._data, UnitSphericalRepresentation)
    else:
        # Check that we have a 2D wrap180 representation
        assert isinstance(hpc1._data, UnitSphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc1.Tx == 0 * u.arcsec
    assert hpc1.Ty == 0 * u.arcsec

    # Check the attrs are in the correct default units
    assert hpc1.Tx.unit is u.arcsec
    assert hpc1.Ty.unit is u.arcsec


@pytest.mark.parametrize(
    'args, kwargs',
    three_D_parameters + [(None, {'Tx': 0 * u.deg,
                                  'Ty': 0 * u.arcsec,
                                  'distance': 1 * u.Mm}),
                          ([0 * u.deg, 0 * u.arcsec], {'distance': 1 * u.Mm})])
def test_create_3d(args, kwargs):
    hpc1 = init_frame(Helioprojective, args, kwargs)

    # Check we have the right class!
    assert isinstance(hpc1, Helioprojective)
    rep_kwarg = kwargs.get('representation', None) if kwargs else None

    if rep_kwarg and rep_kwarg == 'spherical':
        # Check that we have a unitspherical representation
        assert isinstance(hpc1._data, SphericalRepresentation)
    else:
        # Check that we have a 2D wrap180 representation
        assert isinstance(hpc1._data, SphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc1.Tx == 0 * u.arcsec
    assert hpc1.Ty == 0 * u.arcsec
    assert hpc1.distance == 1 * u.Mm

    # Check the attrs are in the correct default units
    assert hpc1.Tx.unit is u.arcsec
    assert hpc1.Ty.unit is u.arcsec
    assert hpc1.distance.unit is u.km


def test_cart_init():
    hpc1 = Helioprojective(CartesianRepresentation(0 * u.km, 0 * u.km, 1 *
                                                   u.Mm))

    assert isinstance(hpc1, Helioprojective)
    assert isinstance(hpc1._data, CartesianRepresentation)


## This is not actually valid, it should be re-purposed to be Heliocentric
## cylindrical Assuming that that is a vaild representation.
#cylindrical_parameters = [
#    ([100 * u.km, 25 * u.deg, 1 * u.Mm], {'representation': 'cylindrical'}), (
#        [100 * u.km, 25 * u.deg, 1 * u.Mm], {'dateobs': '2011/01/01T00:00:00',
#                                             'representation': 'cylindrical'}),
#    ([100 * u.km, 25 * u.deg], {'distance': 1 * u.Mm,
#                                'representation': 'cylindrical'}),
#    (None, {'rho': 100 * u.km,
#            'psi': 25 * u.deg,
#            'distance': 1 * u.Mm,
#            'representation': 'cylindrical'}),
#    ([CylindricalRepresentation(100 * u.km, 25 * u.deg, 1 * u.Mm)],
#     {'representation': 'cylindrical'}), (
#         [CylindricalRepresentation(100 * u.km, 25 * u.deg, 1 * u.Mm)],
#         {'dateobs': '2011/01/01T00:00:00',
#          'representation': 'cylindrical'})
#]
#
#
#@pytest.mark.parametrize('args, kwargs', cylindrical_parameters)
#def test_create_cylindrical(args, kwargs):
#    hpc1 = init_frame(Helioprojective, args, kwargs)
#
#    # Check we have the right class!
#    assert isinstance(hpc1, Helioprojective)
#    # Check that we have a 2D wrap180 representation
#    assert isinstance(hpc1._data, CylindricalRepresentation)
#
#    # Check the attrs are correct
#    assert hpc1.rho == 100 * u.km
#    assert hpc1.psi == 25 * u.deg
#    assert hpc1.distance == 1 * u.Mm
#
#    # Check the attrs are in the correct default units
#    assert hpc1.rho.unit is u.km
#    assert hpc1.psi.unit is u.arcsec
#    assert hpc1.distance.unit is u.km


# Test HPC Calculate Distance
def test_hpc_distance():
    hpc1 = Helioprojective(0 * u.deg, 0 * u.arcsec)

    assert isinstance(hpc1, Helioprojective)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hpc1._data, UnitSphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc1.Tx == 0 * u.arcsec
    assert hpc1.Ty == 0 * u.arcsec

    hpc2 = hpc1.calculate_distance()

    assert isinstance(hpc2._data, SphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc2.Tx == 0 * u.arcsec
    assert hpc2.Ty == 0 * u.arcsec
    assert_quantity_allclose(hpc2.distance, DSUN_METERS - RSUN_METERS)


def test_hpc_distance_off_limb():
    hpc1 = Helioprojective(1500 * u.arcsec, 0 * u.arcsec)

    assert isinstance(hpc1, Helioprojective)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hpc1._data, UnitSphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc1.Tx == 1500 * u.arcsec
    assert hpc1.Ty == 0 * u.arcsec

    hpc2 = hpc1.calculate_distance()

    assert isinstance(hpc2._data, SphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc2.Tx == 1500 * u.arcsec
    assert hpc2.Ty == 0 * u.arcsec
    assert_quantity_allclose(hpc2.distance, u.Quantity(np.nan, u.km))


def test_hpc_distance_3D():
    hpc1 = Helioprojective(1500 * u.arcsec, 0 * u.arcsec, 100 * u.Mm)

    assert isinstance(hpc1, Helioprojective)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hpc1._data, SphericalWrap180Representation)

    # Check the attrs are correct
    assert hpc1.Tx == 1500 * u.arcsec
    assert hpc1.Ty == 0 * u.arcsec

    hpc2 = hpc1.calculate_distance()

    assert hpc2 is hpc1

#==============================================================================
### Heliographic Tests
#==============================================================================


@pytest.mark.parametrize('frame',
                         [HeliographicStonyhurst, HeliographicCarrington])
@pytest.mark.parametrize("args, kwargs", two_D_parameters[:2] + two_D_parameters[4:] + \
                                         [(None, {'lat':0*u.deg, 'lon':0*u.arcsec})])
def test_create_hgs_2d(frame, args, kwargs):
    hgs1 = init_frame(frame, args, kwargs)

    # Check we have the right class!
    assert isinstance(hgs1, frame)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hgs1._data, SphericalWrap180Representation)

    # Check the attrs are correct
    assert hgs1.lon == 0 * u.deg
    assert hgs1.lat == 0 * u.deg

    # Check the attrs are in the correct default units
    assert hgs1.lon.unit is u.deg
    assert hgs1.lat.unit is u.deg
    assert hgs1.radius.unit is u.km


@pytest.mark.parametrize('frame',
                         [HeliographicStonyhurst, HeliographicCarrington])
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
    assert hgs1.lon == 0 * u.deg
    assert hgs1.lat == 0 * u.deg

    # Check the attrs are in the correct default units
    assert hgs1.lon.unit is u.deg
    assert hgs1.lat.unit is u.deg


@pytest.mark.parametrize('frame',
                         [HeliographicStonyhurst, HeliographicCarrington])
@pytest.mark.parametrize(
    "args, kwargs",
    three_D_parameters + [(None, {'lat': 0 * u.deg,
                                  'lon': 0 * u.arcsec,
                                  'radius': 1 * u.Mm}),
                          ([0 * u.deg, 0 * u.arcsec], {'radius': 1 * u.Mm})])
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
    assert hgs1.lon == 0 * u.deg
    assert hgs1.lat == 0 * u.deg
    assert hgs1.radius == 1 * u.Mm

    # Check the attrs are in the correct default units
    assert hgs1.lon.unit is u.deg
    assert hgs1.lat.unit is u.deg
    assert hgs1.radius.unit is u.Mm


def test_hgs_cart_init():
    hpc1 = HeliographicStonyhurst(CartesianRepresentation(0 * u.km, 0 * u.km, 1
                                                          * u.Mm))

    assert isinstance(hpc1, HeliographicStonyhurst)
    assert isinstance(hpc1._data, CartesianRepresentation)


#==============================================================================
### Heliocentric Tests
#==============================================================================
@pytest.mark.parametrize(
    'args, kwargs',
    [((10 * u.km, 10 * u.km, 10 * u.km), None), (None, {'x': 10 * u.km,
                                                        'y': 10 * u.km,
                                                        'z': 10 * u.km}),
     ([CartesianRepresentation(10 * u.km, 10 * u.km, 10 * u.km)], None),
     ([CartesianRepresentation(10 * u.km, 10 * u.km, 10 * u.km)],
      {'dateobs': '2011/01/01T00:00:00'})])
def test_create_hcc_3d(args, kwargs):
    hcc = init_frame(Heliocentric, args, kwargs)

    assert isinstance(hcc, Heliocentric)

    assert isinstance(hcc._data, CartesianRepresentation)

    assert hcc.x == 10 * u.km
    assert hcc.y == 10 * u.km
    assert hcc.z == 10 * u.km

    # Check the attrs are in the correct default units
    assert hcc.x.unit is u.km
    assert hcc.y.unit is u.km
    assert hcc.z.unit is u.km
