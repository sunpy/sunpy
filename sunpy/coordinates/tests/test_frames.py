
import warnings
from contextlib import nullcontext

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
from sunpy.coordinates import PlanarScreen, SphericalScreen, propagate_with_solar_surface
from sunpy.coordinates.frames import (
    Geomagnetic,
    Heliocentric,
    HeliographicCarrington,
    HeliographicStonyhurst,
    Helioprojective,
    HelioprojectiveRadial,
)
from sunpy.coordinates.sun import angular_radius
from sunpy.time import parse_time
from sunpy.util.exceptions import SunpyDeprecationWarning, SunpyUserWarning

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


@pytest.mark.parametrize(('args', 'kwargs'),
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
    ('args', 'kwargs'),
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

    with pytest.warns(SunpyUserWarning, match="is all NaNs"):
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
    with np.errstate(over='ignore'):
        hpc = Helioprojective(u.Quantity(0, u.deg, dtype=np.float32),
                              u.Quantity(0, u.arcsec, dtype=np.float16),
                              observer=HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU))

        with pytest.warns(SunpyUserWarning, match="Tx is float32, and Ty is float16"):
            hpc.make_3d()


def test_hpc_obstime_from_observer():
    # Test that observer.obstime is used for obstime if obstime is not provided
    observer = HeliographicStonyhurst(obstime='2023-09-08')
    hpc = Helioprojective(observer=observer)
    assert hpc.obstime == observer.obstime

    # Test that obstime is None if observer does not have an obstime
    hpc = Helioprojective(observer='earth')
    assert hpc.obstime is None


def test_hpc_is_visible_2d():
    hpc = Helioprojective(2000*u.arcsec, 2000*u.arcsec,
                          observer='earth', obstime='2023-08-03')
    assert hpc.is_visible()


def test_hpc_is_visible():
    hpc = Helioprojective([0]*2*u.arcsec, [0]*2*u.arcsec, [0.5, 1.5]*u.AU,
                          observer='earth', obstime='2023-08-03')
    assert (hpc.is_visible() == [True, False]).all()


def test_hpc_is_visible_tolerance():
    hpc = Helioprojective(200*u.arcsec, 0*u.arcsec,
                          observer='earth', obstime='2023-08-03').make_3d()

    # Due to the limitations of numerical precision, the coordinate may be computed to be slightly
    # below the solar surface, and thus may be invisible when the tolerance is set to zero
    if hpc.is_visible(tolerance=0*u.m):
        pytest.skip("Test already passes prior to increasing the tolerance.")

    assert hpc.is_visible(tolerance=1*u.m)


# ==============================================================================
# Helioprojective Radial Tests
# ==============================================================================

@pytest.mark.parametrize(('args', 'kwargs'),
                         two_D_parameters + [(None, {'psi': 0 * u.deg,
                                                     'delta': 0 * u.arcsec})])
def test_create_hpr_2d(args, kwargs):
    hpr1 = init_frame(HelioprojectiveRadial, args, kwargs)

    assert isinstance(hpr1, HelioprojectiveRadial)
    assert isinstance(hpr1._data, UnitSphericalRepresentation)

    assert hpr1.psi.unit is u.deg
    assert hpr1.delta.unit is u.deg
    assert_quantity_allclose(hpr1.psi, 0*u.deg)
    assert_quantity_allclose(hpr1.delta, 0*u.deg)

    assert_quantity_allclose(hpr1.theta, 90*u.deg)


@pytest.mark.parametrize(
    ('args', 'kwargs'),
    three_D_parameters + [(None, {'psi': 0 * u.deg,
                                  'delta': 0 * u.arcsec,
                                  'distance': 1 * u.Mm}),
                          ([0 * u.deg, 0 * u.arcsec], {'distance': 1 * u.Mm})])
def test_create_hpr_3d(args, kwargs):
    hpr1 = init_frame(HelioprojectiveRadial, args, kwargs)

    assert isinstance(hpr1, HelioprojectiveRadial)
    assert isinstance(hpr1._data, SphericalRepresentation)

    assert hpr1.psi.unit is u.deg
    assert hpr1.delta.unit is u.deg
    assert hpr1.distance.unit is u.Mm
    assert_quantity_allclose(hpr1.psi, 0*u.deg)
    assert_quantity_allclose(hpr1.delta, 0*u.deg)
    assert_quantity_allclose(hpr1.distance, 1*u.Mm)

    assert_quantity_allclose(hpr1.theta, 90*u.deg)

    # Since hpr1 is already 3D, make_3d() should simply return the original object
    hpr2 = hpr1.make_3d()
    assert hpr2 is hpr1


def test_hpr_distance():
    hpr1 = HelioprojectiveRadial(0*u.deg, -90*u.deg,
                                 observer=HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU))

    hpr2 = hpr1.make_3d()

    assert_quantity_allclose(hpr2.psi, 0*u.deg)
    assert_quantity_allclose(hpr2.delta, -90*u.deg)
    assert_quantity_allclose(hpr2.distance, DSUN_METERS - RSUN_METERS)


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
@pytest.mark.parametrize(('args', 'kwargs'), two_D_parameters +
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
    ('args', 'kwargs'),
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
    ('args', 'kwargs'),
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


@pytest.mark.parametrize(('x', 'y', 'psi'), [(0*u.km, -1*u.km, 270*u.deg),
                                             (0*u.km, 1*u.km, 90*u.deg),
                                             (-1*u.km, 0*u.km, 180*u.deg)])
def test_heliocentric_radial_psi(x, y, psi):
    # The cylindrical representation of HCC is Heliocentric Radial
    # Test that the `psi` component is represented as desired
    # The definition is shifted by 90 degrees relative to Thompson (2006)
    hcc = Heliocentric(CartesianRepresentation(x, y, 0*u.km), representation_type='cylindrical')

    assert_quantity_allclose(hcc.psi, psi)


# ==============================================================================
# Magnetic-model coordinate frame tests
# ==============================================================================


@pytest.mark.parametrize(('obstime', 'dipole_lonlat'),
                         [('2012-07-01', [-72.408328, 80.16423]*u.deg),
                          ('2012-08-01', [-72.415148, 80.169261]*u.deg),
                          ('2032-08-01', [-72.576963, 81.212062]*u.deg),
                          (['2012-07-01', '2012-08-01'], [[-72.408328, 80.16423], [-72.415148, 80.169261]]*u.deg),
                          (['2012-07-01', '2032-08-01'], [[-72.408328, 80.16423], [-72.576963, 81.212062]]*u.deg),
                          ('1899-01-01', ValueError),
                          (['1899-01-01', '2012-07-01'], ValueError)])
def test_magnetic_model_obstime(obstime, dipole_lonlat):
    frame = Geomagnetic(obstime=obstime, magnetic_model='igrf13')
    ctx = pytest.raises(ValueError, match="earlier than the year 1900") if dipole_lonlat is ValueError else nullcontext()
    with ctx:
        assert_quantity_allclose(frame.dipole_lonlat, dipole_lonlat.T)


def test_magnetic_model_default():
    # Also tests that no downloading happens because this test is not marked as remote
    obstime = '2012-07-01'
    frame_default = Geomagnetic(obstime=obstime)
    frame_igrf13 = Geomagnetic(obstime=obstime, magnetic_model='igrf13')

    assert_quantity_allclose(frame_igrf13.dipole_lonlat, [-72.408328, 80.16423]*u.deg)
    assert_quantity_allclose(frame_default.dipole_lonlat, frame_igrf13.dipole_lonlat)


@pytest.mark.remote_data
@pytest.mark.parametrize(('magnetic_model', 'obstime', 'dipole_lonlat'),
                         [('igrf12', '2012-07-01', [-72.414318, 80.16354]*u.deg),
                          ('igrf11', '2006-01-01', [-71.886023, 79.801523]*u.deg),
                          ('igrf10', '2006-01-01', [-71.822653, 79.785185]*u.deg)])
def test_magnetic_model_with(magnetic_model, obstime, dipole_lonlat):
    frame = Geomagnetic(magnetic_model=magnetic_model, obstime=obstime)
    assert_quantity_allclose(frame.dipole_lonlat, dipole_lonlat)


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


@pytest.mark.parametrize(('args', 'kwargs'),
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
    # Only 1 of the 7 parameterized arguments actually emits the warning
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=(
            'The conversion of these 2D helioprojective coordinates to '
        ), category=SunpyUserWarning)
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


@pytest.fixture
def off_limb_coord():
    frame = Helioprojective(observer='earth', obstime='2020-01-01')
    return SkyCoord(Tx=[-1000, 300, 1000]*u.arcsec, Ty=[-1000, 300, 1000]*u.arcsec, frame=frame)


@pytest.fixture
def non_earth_coord():
    return SkyCoord(70*u.deg, 20*u.deg, 1*u.AU, obstime='2020-01-01', frame=HeliographicStonyhurst)


@pytest.mark.parametrize('screen_class', [
    SphericalScreen,
    PlanarScreen,
])
def test_screen_classes(off_limb_coord, screen_class):
    # Smoke test for spherical screen
    with pytest.warns(SunpyUserWarning, match='The conversion of these 2D helioprojective coordinates to 3D is all NaNs'):
            olc_3d = off_limb_coord[0].make_3d()
    assert np.isnan(olc_3d.distance).all()
    sph_screen = screen_class(off_limb_coord[0].observer)
    with sph_screen:
        olc_3d = off_limb_coord[0].make_3d()
    assert not np.isnan(olc_3d.distance).all()


def test_assume_spherical_screen_deprecated(off_limb_coord):
    with pytest.warns(SunpyDeprecationWarning, match='The assume_spherical_screen function is deprecated'):
        with Helioprojective.assume_spherical_screen(off_limb_coord.observer):
            _ = off_limb_coord.make_3d()


@pytest.mark.parametrize(('only_off_disk', 'distance_from_center', 'distance'), [
    (False, 0*u.m, [0.98331616, 0.98329512, 0.98331616]*u.AU),
    (True, 0*u.m, [0.98331616, 0.97910333, 0.98331616]*u.AU),
    (False, 1*u.Rsun, [0.97866558, 0.97864465, 0.97866558]*u.AU),
])
def test_planar_screen(off_limb_coord, only_off_disk, distance_from_center, distance):
    with PlanarScreen(off_limb_coord.observer, distance_from_center=distance_from_center, only_off_disk=only_off_disk):
        olc_3d = off_limb_coord.make_3d()
    assert u.quantity.allclose(olc_3d.distance, distance)


@pytest.mark.parametrize(('only_off_disk', 'distance'), [
    (False, [0.98329304, 0.98329304, 0.98329304]*u.AU),
    (True, [0.98329304, 0.97910333, 0.98329304]*u.AU),
])
def test_spherical_screen(off_limb_coord, only_off_disk, distance):
    with SphericalScreen(off_limb_coord.observer, only_off_disk=only_off_disk):
        olc_3d = off_limb_coord.make_3d()
    assert u.quantity.allclose(olc_3d.distance, distance)


@pytest.mark.parametrize(('only_off_disk', 'distance_from_center', 'distance'), [
    (False, 0*u.m, [0.96419505, 0.98918004, 1.00321100]*u.AU),
    (True, 0*u.m, [0.96419505, 0.97910333, 1.00321100]*u.AU),
    (False, 1*u.Rsun, [0.94916557, 0.97376111, 0.98757335]*u.AU),
])
def test_planar_screen_askew(off_limb_coord, only_off_disk, distance_from_center, distance, non_earth_coord):
    with PlanarScreen(non_earth_coord, distance_from_center=distance_from_center, only_off_disk=only_off_disk):
        olc_3d = off_limb_coord.make_3d()
    assert u.quantity.allclose(olc_3d.distance, distance)


@pytest.mark.parametrize(('only_off_disk', 'distance'), [
    (False, [0.96348934, 0.98911699, 1.00251206]*u.AU),
    (True, [0.96348934, 0.97910333, 1.00251206]*u.AU),
])
def test_spherical_screen_askew(off_limb_coord, only_off_disk, distance, non_earth_coord):
    with SphericalScreen(non_earth_coord, only_off_disk=only_off_disk):
        olc_3d = off_limb_coord.make_3d()
    assert u.quantity.allclose(olc_3d.distance, distance)


@pytest.mark.parametrize(('screen', 'only_off_disk', 'distance'), [
    (SphericalScreen, False, [0.98405002, 0.98306592, 0.98253594]*u.AU),
    (SphericalScreen, True, [0.98405002, 0.97910333, 0.98253594]*u.AU),
    (PlanarScreen, False, [0.98407381, 0.98306806, 0.98255967]*u.AU),
    (PlanarScreen, True, [0.98407381, 0.97910333, 0.98255967]*u.AU),
])
def test_screen_plus_diffrot(off_limb_coord, screen, only_off_disk, distance):
    new_observer = HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU, obstime=off_limb_coord.obstime + 1*u.day)
    with propagate_with_solar_surface(), screen(new_observer, only_off_disk=only_off_disk):
        olc_3d = off_limb_coord.make_3d()
    assert_quantity_allclose(olc_3d.distance, distance)


@pytest.mark.parametrize('only_off_disk', [False, True])
@pytest.mark.parametrize('screen', [SphericalScreen, PlanarScreen])
def test_screen_plus_diffrot_array_obstime(screen, only_off_disk):
    array_obstime = parse_time('2025-05-30') + np.arange(5) * u.day
    observer = HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU, obstime=array_obstime)
    coord_2d = SkyCoord([-1000, -750, -500, -250, 0] * u.arcsec,
                        [600]*5 * u.arcsec,
                        frame=Helioprojective(observer=observer, obstime=array_obstime))

    with propagate_with_solar_surface(), screen(observer[0], only_off_disk=only_off_disk):
        # Confirm that array obstime works
        coord_3d = coord_2d.make_3d()

        # Confirm that calculated distances match the calculation when done individually
        for point_2d, point_3d in zip(coord_2d, coord_3d):
            point_2d_3d = point_2d.make_3d()
            assert_quantity_allclose(point_2d_3d.separation_3d(point_3d), 0*u.m, atol=1*u.m)
