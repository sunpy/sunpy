
import warnings

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import (
    CartesianRepresentation,
    SkyCoord,
    SphericalRepresentation,
    UnitSphericalRepresentation,
)
from astropy.tests.helper import assert_quantity_allclose

from sunpy import sun
from sunpy.coordinates.frames import (
    Heliocentric,
    HeliographicCarrington,
    HeliographicStonyhurst,
    Helioprojective,
)
from sunpy.coordinates.sun import angular_radius
from sunpy.time import parse_time
from sunpy.util.exceptions import SunpyUserWarning

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
    ([0 * u.deg, 0 * u.arcsec], {'obstime': '2011/01/01T00:00:00'}),
    ([UnitSphericalRepresentation(0 * u.deg, 0 * u.arcsec)], None),
    ([UnitSphericalRepresentation(0 * u.deg, 0 * u.arcsec)], {'obstime': '2011/01/01T00:00:00'}),
    ([0 * u.deg, 0 * u.arcsec], {'representation_type': 'unitspherical'})
]
"""
These are common 3D params, kwargs are frame specific
"""
three_D_parameters = [
    ([0 * u.deg, 0 * u.arcsec, 1 * u.Mm], None),
    ([0 * u.deg, 0 * u.arcsec, 1 * u.Mm], {'obstime': '2011/01/01T00:00:00'}),
    ([0 * u.deg, 0 * u.arcsec, 1 * u.Mm], {'representation_type': 'spherical'}),
    ([SphericalRepresentation(0 * u.deg, 0 * u.arcsec, 1 * u.Mm)],
     None),
    ([SphericalRepresentation(0 * u.deg, 0 * u.arcsec, 1 * u.Mm)], None), (
        [SphericalRepresentation(0 * u.deg, 0 * u.arcsec, 1 * u.Mm)],
        {'obstime': '2011/01/01T00:00:00'})
]

# ==============================================================================
# Helioprojective Tests
# ==============================================================================


@pytest.mark.parametrize('args, kwargs',
                         two_D_parameters + [(None, {'Tx': 0 * u.deg,
                                                     'Ty': 0 * u.arcsec})])
def test_create_hpc_2d(args, kwargs):
    hpc1 = init_frame(Helioprojective, args, kwargs)

    # Check we have the right class!
    assert isinstance(hpc1, Helioprojective)
    rep_kwarg = kwargs.get('representation_type', None) if kwargs else None

    if rep_kwarg and rep_kwarg == 'unitspherical':
        # Check that we have a unitspherical representation
        assert isinstance(hpc1._data, UnitSphericalRepresentation)
    else:
        # Check that we have a 2D wrap180 representation
        assert isinstance(hpc1._data, UnitSphericalRepresentation)

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
    rep_kwarg = kwargs.get('representation_type', None) if kwargs else None

    if rep_kwarg and rep_kwarg == 'spherical':
        # Check that we have a unitspherical representation
        assert isinstance(hpc1._data, SphericalRepresentation)
    else:
        # Check that we have a 2D wrap180 representation
        assert isinstance(hpc1._data, SphericalRepresentation)

    # Check the attrs are correct
    assert hpc1.Tx == 0 * u.arcsec
    assert hpc1.Ty == 0 * u.arcsec
    assert hpc1.distance == 1 * u.Mm

    # Check the attrs are in the correct default units
    assert hpc1.Tx.unit is u.arcsec
    assert hpc1.Ty.unit is u.arcsec
    assert hpc1.distance.unit is u.Mm


def test_cart_init():
    hpc1 = Helioprojective(CartesianRepresentation(0 * u.km, 0 * u.km, 1 *
                                                   u.Mm))

    assert isinstance(hpc1, Helioprojective)
    assert isinstance(hpc1._data, CartesianRepresentation)


# Test HPC Calculate Distance
def test_hpc_distance():
    hpc1 = Helioprojective(0 * u.deg, 0 * u.arcsec,
                           observer=HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU))

    assert isinstance(hpc1, Helioprojective)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hpc1._data, UnitSphericalRepresentation)

    # Check the attrs are correct
    assert hpc1.Tx == 0 * u.arcsec
    assert hpc1.Ty == 0 * u.arcsec

    hpc2 = hpc1.make_3d()

    assert isinstance(hpc2._data, SphericalRepresentation)

    # Check the attrs are correct
    assert hpc2.Tx == 0 * u.arcsec
    assert hpc2.Ty == 0 * u.arcsec
    assert_quantity_allclose(hpc2.distance, DSUN_METERS - RSUN_METERS)


def test_hpc_distance_cartesian():
    # Test detection of distance in other representations
    hpc1 = Helioprojective(CartesianRepresentation(0 * u.km, 0 * u.km, 1 * u.Mm))

    assert isinstance(hpc1, Helioprojective)
    assert isinstance(hpc1._data, CartesianRepresentation)

    assert hpc1.make_3d() is hpc1


def test_hpc_distance_off_limb():
    hpc1 = Helioprojective(1500 * u.arcsec, 0 * u.arcsec,
                           observer=HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU))

    assert isinstance(hpc1, Helioprojective)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hpc1._data, UnitSphericalRepresentation)

    # Check the attrs are correct
    assert hpc1.Tx == 1500 * u.arcsec
    assert hpc1.Ty == 0 * u.arcsec

    hpc2 = hpc1.make_3d()

    assert isinstance(hpc2._data, SphericalRepresentation)

    # Check the attrs are correct
    assert hpc2.Tx == 1500 * u.arcsec
    assert hpc2.Ty == 0 * u.arcsec
    assert_quantity_allclose(hpc2.distance, u.Quantity(np.nan, u.km))


def test_hpc_distance_3D():
    hpc1 = Helioprojective(1500 * u.arcsec, 0 * u.arcsec, 100 * u.Mm)

    assert isinstance(hpc1, Helioprojective)
    # Check that we have a 2D wrap180 representation
    assert isinstance(hpc1._data, SphericalRepresentation)

    # Check the attrs are correct
    assert hpc1.Tx == 1500 * u.arcsec
    assert hpc1.Ty == 0 * u.arcsec

    hpc2 = hpc1.make_3d()

    assert hpc2 is hpc1


def test_wrapping_on():
    hpc1 = Helioprojective(359.9*u.deg, 10*u.deg)
    assert_quantity_allclose(hpc1.Tx, -0.1*u.deg)
    assert_quantity_allclose(hpc1.Tx.wrap_angle, 180*u.deg)


def test_wrapping_off():
    hpc1 = Helioprojective(359.9*u.deg, 10*u.deg, wrap_longitude=False)
    assert_quantity_allclose(hpc1.Tx, 359.9*u.deg)
    assert_quantity_allclose(hpc1.Tx.wrap_angle, 360*u.deg)


def test_hpc_default_observer():
    # Observer is considered default if it hasn't been specified *and* if obstime isn't specified
    hpc = Helioprojective(0*u.arcsec, 0*u.arcsec)
    assert hpc.is_frame_attr_default('observer')

    hpc = Helioprojective(0*u.arcsec, 0*u.arcsec, obstime='2019-06-01')
    assert not hpc.is_frame_attr_default('observer')


def test_hpc_low_precision_float_warning():
    hpc = Helioprojective(u.Quantity(0, u.deg, dtype=np.float32),
                          u.Quantity(0, u.arcsec, dtype=np.float16),
                          observer=HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU))

    with pytest.raises(SunpyUserWarning, match="Tx is float32, and Ty is float16"):
        hpc.make_3d()


# ==============================================================================
# ## Heliographic Tests
# ==============================================================================

def test_HEEQ_creation():
    # Smoke test to make sure HEEQ constructors work fine
    _ = HeliographicStonyhurst(lon=0*u.deg, lat=90*u.deg,
                               obstime=parse_time('2018-12-21'))
    _ = HeliographicStonyhurst(lon=0*u.deg, lat=90*u.deg, radius=1*u.km,
                               obstime=parse_time('2018-12-21'))
    _ = HeliographicStonyhurst(x=1*u.km, y=1*u.km, z=1*u.km,
                               obstime=parse_time('2018-12-21'),
                               representation_type='cartesian')


@pytest.mark.parametrize('frame',
                         [HeliographicStonyhurst, HeliographicCarrington])
@pytest.mark.parametrize("args, kwargs", two_D_parameters +
                         [(None, {'lat': 0*u.deg, 'lon': 0*u.arcsec})])
def test_create_hgs_2d(frame, args, kwargs):
    hgs1 = init_frame(frame, args, kwargs)

    # Check we have the right class!
    assert isinstance(hgs1, frame)
    # Check that we have a 2D representation
    assert isinstance(hgs1._data, UnitSphericalRepresentation)

    # Check the attrs are correct
    assert hgs1.lon == 0 * u.deg
    assert hgs1.lat == 0 * u.deg

    # Check the attrs are in the correct default units
    assert hgs1.lon.unit is u.deg
    assert hgs1.lat.unit is u.deg

    # Test the value of the rsun frame attribute
    assert_quantity_allclose(hgs1.rsun, sun.constants.radius)

    # Test conversion to 3D
    hgs_3d = hgs1.make_3d()
    assert_quantity_allclose(hgs_3d.lon, hgs1.lon)
    assert_quantity_allclose(hgs_3d.lat, hgs1.lat)
    assert_quantity_allclose(hgs_3d.radius, hgs1.rsun)


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

    rep_kwarg = kwargs.get('representation_type', None) if kwargs else None

    if rep_kwarg == 'spherical':
        assert isinstance(hgs1._data, SphericalRepresentation)
    else:
        # Check Carrington first because it's a subclass of Stonyhurst
        if isinstance(hgs1, HeliographicCarrington):
            # Check that we have a 2D wrap180 representation
            assert isinstance(hgs1._data, SphericalRepresentation)
        elif isinstance(hgs1, HeliographicStonyhurst):
            # Check that we have a 2D wrap180 representation
            assert isinstance(hgs1._data, SphericalRepresentation)

    # Check the attrs are correct
    assert hgs1.lon == 0 * u.deg
    assert hgs1.lat == 0 * u.deg
    assert hgs1.radius == 1 * u.Mm

    # Check the attrs are in the correct default units
    assert hgs1.lon.unit is u.deg
    assert hgs1.lat.unit is u.deg
    assert hgs1.radius.unit is u.Mm


def test_hgs_cart_init():
    hpc1 = HeliographicStonyhurst(CartesianRepresentation(0 * u.km,
                                                          0 * u.km,
                                                          1 * u.Mm))

    assert isinstance(hpc1, HeliographicStonyhurst)
    assert isinstance(hpc1._data, CartesianRepresentation)


def test_hgs_wrapping_on():
    hpc1 = HeliographicStonyhurst(350*u.deg, 10*u.deg)
    assert_quantity_allclose(hpc1.lon, -10*u.deg)
    assert_quantity_allclose(hpc1.lon.wrap_angle, 180*u.deg)


def test_hgs_wrapping_off():
    hpc1 = HeliographicStonyhurst(350*u.deg, 10*u.deg, wrap_longitude=False)
    assert_quantity_allclose(hpc1.lon, 350*u.deg)
    assert_quantity_allclose(hpc1.lon.wrap_angle, 360*u.deg)


def test_hgc_wrapping_360():
    hpc1 = HeliographicCarrington(350*u.deg, 10*u.deg)
    assert_quantity_allclose(hpc1.lon, 350*u.deg)
    assert_quantity_allclose(hpc1.lon.wrap_angle, 360*u.deg)


# ==============================================================================
# ## Heliocentric Tests
# ==============================================================================
@pytest.mark.parametrize(
    'args, kwargs',
    [((10 * u.km, 10 * u.km, 10 * u.km), None), (None, {'x': 10 * u.km,
                                                        'y': 10 * u.km,
                                                        'z': 10 * u.km}),
     ([CartesianRepresentation(10 * u.km, 10 * u.km, 10 * u.km)], None),
     ([CartesianRepresentation(10 * u.km, 10 * u.km, 10 * u.km)],
      {'obstime': '2011/01/01T00:00:00'})])
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


def test_hcc_default_observer():
    # Observer is considered default if it hasn't been specified *and* if obstime isn't specified
    hcc = Heliocentric(0*u.AU, 0*u.AU, 0*u.AU)
    assert hcc.is_frame_attr_default('observer')

    hcc = Heliocentric(0*u.AU, 0*u.AU, 0*u.AU, obstime='2019-06-01')
    assert not hcc.is_frame_attr_default('observer')


@pytest.mark.parametrize('x, y, psi', [(0*u.km, -1*u.km, 270*u.deg),
                                       (0*u.km, 1*u.km, 90*u.deg),
                                       (-1*u.km, 0*u.km, 180*u.deg)])
def test_heliocentric_radial_psi(x, y, psi):
    # The cylindrical representation of HCC is Heliocentric Radial
    # Test that the `psi` component is represented as desired
    # The definition is shifted by 90 degrees relative to Thompson (2006)
    hcc = Heliocentric(CartesianRepresentation(x, y, 0*u.km), representation_type='cylindrical')

    assert_quantity_allclose(hcc.psi, psi)


# ==============================================================================
# SkyCoord Tests
# ==============================================================================


two_D_parameters = [
    ([0 * u.deg, 0 * u.arcsec], {}),
    ([[0, 1, 2, 3], [5, 6, 7, 8]], {'unit': u.deg}),
    ([0 * u.deg, 0 * u.arcsec], {'representation_type': SphericalRepresentation}),
    ([UnitSphericalRepresentation(0 * u.deg, 0 * u.arcsec)], {}),
    ([UnitSphericalRepresentation(0 * u.deg, 0 * u.arcsec)], {}),
    ([SphericalRepresentation(0 * u.deg, 0 * u.arcsec, 1*u.one)], {}),
]


@pytest.mark.parametrize("args, kwargs",
                         two_D_parameters + [([0 * u.deg, 0 * u.arcsec],
                                              {'representation_type': 'unitspherical'})])
def test_skycoord_hpc(args, kwargs):
    """
    Test that when instantiating a HPC frame with SkyCoord that make_3d
    still works.
    """

    sc = SkyCoord(*args, **kwargs, frame="helioprojective",
                  observer='earth', obstime="2011-01-01T00:00:00")
    # Test the transform to HGS because it will force a `make_3d` call.
    #
    # astropy emits some warnings here because of invalid NaN comparisons,
    # which will be removed in a future astropy release
    # (see https://github.com/astropy/astropy/pull/9843), so filter the warnings out.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='invalid value encountered',
                                category=RuntimeWarning)
        hgs = sc.transform_to("heliographic_stonyhurst")

    assert isinstance(hgs.frame, HeliographicStonyhurst)


def test_hgc_incomplete_observer():
    with pytest.raises(ValueError, match=r'Full 3D coordinate \(including radius\) must be specified'):
        SkyCoord(0*u.deg, 0*u.deg, frame="heliographic_carrington",
                 observer='self', obstime="2011-01-01T00:00:00")


def test_angular_radius():
    coord = Helioprojective(0*u.deg, 0*u.deg, 5*u.km, obstime="2010/01/01T00:00:00", observer="earth")
    assert_quantity_allclose(coord.angular_radius, angular_radius(coord.obstime))


def test_angular_radius_no_observer():
    coord = Helioprojective(0*u.deg, 0*u.deg, 5*u.km, obstime="2010/01/01T00:00:00", observer=None)
    with pytest.raises(ValueError, match=r"The observer must be defined, not `None`"):
        coord.angular_radius


def test_angular_radius_no_obstime():
    coord = Helioprojective(0*u.deg, 0*u.deg, 5*u.km, obstime=None, observer="earth")
    with pytest.raises(ValueError, match=r"The observer must be fully defined by specifying `obstime`."):
        coord.angular_radius
