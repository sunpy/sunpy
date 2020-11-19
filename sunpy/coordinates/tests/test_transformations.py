import numpy as np
import pytest

import astropy.units as u
from astropy.constants import c as speed_of_light
from astropy.coordinates import (
    ICRS,
    Angle,
    CartesianDifferential,
    CartesianRepresentation,
    ConvertError,
    HeliocentricMeanEcliptic,
    Longitude,
    SkyCoord,
    SphericalDifferential,
    get_body_barycentric,
    get_body_barycentric_posvel,
)
from astropy.tests.helper import assert_quantity_allclose, quantity_allclose
from astropy.time import Time

from sunpy.coordinates import (
    GeocentricEarthEquatorial,
    GeocentricSolarEcliptic,
    Heliocentric,
    HeliocentricEarthEcliptic,
    HeliocentricInertial,
    HeliographicCarrington,
    HeliographicStonyhurst,
    Helioprojective,
    sun,
)
from sunpy.coordinates.ephemeris import get_body_heliographic_stonyhurst, get_earth
from sunpy.coordinates.frames import _J2000
from sunpy.coordinates.transformations import transform_with_sun_center
from sunpy.sun.constants import radius as _RSUN
from sunpy.sun.constants import sidereal_rotation_rate
from sunpy.time import parse_time


def test_hcc_to_hgs():
    '''
    Check that a coordinate pointing to the observer in Heliocentric
    coordinates maps to the lattitude/longitude of the observer in
    HeliographicStonyhurst coordinates.
    '''
    lat = 10 * u.deg
    lon = 20 * u.deg
    observer = HeliographicStonyhurst(lat=lat, lon=lon)
    hcc_in = Heliocentric(x=0*u.km, y=0*u.km, z=1*u.km, observer=observer)
    hgs_out = hcc_in.transform_to(HeliographicStonyhurst())

    assert_quantity_allclose(hgs_out.lat, lat)
    assert_quantity_allclose(hgs_out.lon, lon)


def test_hpc_hpc():
    # Use some unphysical values for solar parameters for testing, to make it
    # easier to calculate expected results.
    rsun = 1*u.m
    D0 = 1*u.km
    L0 = 1*u.deg
    observer_in = HeliographicStonyhurst(lat=0*u.deg, lon=0*u.deg, radius=D0)
    observer_out = HeliographicStonyhurst(lat=0*u.deg, lon=L0, radius=D0)

    hpc_in = Helioprojective(0*u.arcsec, 0*u.arcsec, rsun=rsun, observer=observer_in)
    hpc_out = Helioprojective(observer=observer_out, rsun=rsun)

    hpc_new = hpc_in.transform_to(hpc_out)

    assert hpc_new.observer == hpc_out.observer

    # Calculate the distance subtended by an angle of L0 from the centre of the
    # Sun.
    dd = -1 * rsun * np.tan(L0)
    # Calculate the angle corresponding to that distance as seen by the new
    # observer.
    theta = np.arctan2(dd, (D0 - rsun))

    assert quantity_allclose(theta, hpc_new.Tx, rtol=1e-3)


def test_hpc_hpc_sc():
    # Use some unphysical values for solar parameters for testing, to make it
    # easier to calculate expected results.
    rsun = 1*u.m
    D0 = 1*u.km
    L0 = 1*u.deg
    observer_in = HeliographicStonyhurst(lat=0*u.deg, lon=0*u.deg, radius=D0)
    observer_out = HeliographicStonyhurst(lat=0*u.deg, lon=L0, radius=D0)

    sc_in = SkyCoord(0*u.arcsec, 0*u.arcsec, rsun=rsun, observer=observer_in,
                     frame='helioprojective')
    hpc_out = Helioprojective(observer=observer_out, rsun=rsun)

    hpc_new = sc_in.transform_to(hpc_out)

    assert hpc_new.observer.lat == hpc_out.observer.lat
    assert hpc_new.observer.lon == hpc_out.observer.lon
    assert hpc_new.observer.radius == hpc_out.observer.radius


def test_hpc_hpc_null():
    hpc_in = Helioprojective(0*u.arcsec, 0*u.arcsec)
    hpc_out = Helioprojective()

    hpc_new = hpc_in.transform_to(hpc_out)

    assert hpc_new is not hpc_in
    assert quantity_allclose(hpc_new.Tx, hpc_in.Tx)
    assert quantity_allclose(hpc_new.Ty, hpc_in.Ty)
    assert hpc_out.observer == hpc_new.observer


def test_hcrs_hgs():
    # Get the current Earth location in HCRS
    adate = parse_time('2015/05/01 01:13:00')
    earth_hcrs = SkyCoord(get_body_barycentric('earth', adate), frame='icrs', obstime=adate).hcrs

    # Convert from HCRS to HGS
    earth_hgs = earth_hcrs.transform_to(HeliographicStonyhurst)

    # The HGS longitude of the Earth should be zero within numerical error
    # Due to an issue with wrapping at +-360, we shift it to pass the test.
    assert quantity_allclose((earth_hgs.lon+1*u.deg) % (360*u.deg), 1*u.deg, atol=1e-12*u.deg)

    # The HGS latitude and radius should be within valid ranges
    assert quantity_allclose(earth_hgs.lat, 0*u.deg, atol=7.3*u.deg)
    assert quantity_allclose(earth_hgs.radius, 1*u.AU, atol=0.017*u.AU)


def test_hcrs_hgs_array_obstime():
    # Get the Earth location in HCRS at two times
    times = Time(['2017-01-01', '2017-06-01'])
    earth_hcrs = SkyCoord(get_body_barycentric('earth', times), frame='icrs', obstime=times).hcrs

    # Transform each time in separate calls (uses scalar obstime)
    earth_hgs_0 = earth_hcrs[0].transform_to(HeliographicStonyhurst)
    earth_hgs_1 = earth_hcrs[1].transform_to(HeliographicStonyhurst)

    # Transform both times in one call (uses array obstime)
    earth_hgs = earth_hcrs.transform_to(HeliographicStonyhurst)

    # Confirm that the two approaches produce the same results
    assert quantity_allclose(earth_hgs_0.lon, earth_hgs[0].lon, atol=1e-12*u.deg)
    assert quantity_allclose(earth_hgs_0.lat, earth_hgs[0].lat, rtol=1e-10)
    assert quantity_allclose(earth_hgs_0.radius, earth_hgs[0].radius, rtol=1e-10)
    assert quantity_allclose(earth_hgs_1.lon, earth_hgs[1].lon, atol=1e-12*u.deg)
    assert quantity_allclose(earth_hgs_1.lat, earth_hgs[1].lat, rtol=1e-10)
    assert quantity_allclose(earth_hgs_1.radius, earth_hgs[1].radius, rtol=1e-10)


def test_hgs_hcrs():
    # This test checks the HGS->HCRS transformation by transforming from HGS to
    # HeliocentricMeanEcliptic (HME).  It will fail if there are errors in Astropy's
    # HCRS->ICRS or ICRS->HME transformations.

    # Use published HGS coordinates in the Astronomical Almanac (2013), pages C6-C7
    obstime = Time('2013-01-28')
    earth_hgs = SkyCoord(0*u.deg, -5.73*u.deg, 0.9848139*u.AU, frame=HeliographicStonyhurst,
                         obstime=obstime)

    # Transform to HME at observation-time equinox
    earth_hme = earth_hgs.transform_to(HeliocentricMeanEcliptic(equinox=obstime))

    # Validate against published values from the Astronomical Almanac (2013), page C6 per page E2
    # The dominant source of inaccuracy is the limited precision of the published B0 used above
    assert quantity_allclose(earth_hme.lon, Angle('308d13m30.51s') - 180*u.deg, atol=5*u.arcsec)
    assert quantity_allclose(earth_hme.lat, -Angle('-0.27s'), atol=10*u.arcsec)
    assert quantity_allclose(earth_hme.distance, 0.9848139*u.AU, atol=5e-7*u.AU)


def test_hgs_hgc_roundtrip():
    obstime = "2011-01-01"

    hgsin = HeliographicStonyhurst(lat=10*u.deg, lon=20*u.deg, obstime=obstime)
    hgcout = hgsin.transform_to(HeliographicCarrington(observer='earth', obstime=obstime))

    assert_quantity_allclose(hgsin.lat, hgcout.lat)
    assert_quantity_allclose(hgsin.lon + sun.L0(obstime), hgcout.lon)

    hgsout = hgcout.transform_to(HeliographicStonyhurst(obstime=obstime))

    assert_quantity_allclose(hgsout.lat, hgsin.lat)
    assert_quantity_allclose(hgsout.lon, hgsin.lon)


def test_hgs_cartesian_rep_to_hpc():
    # This test checks transformation HGS->HPC when the coordinate is in a Cartesian
    # representation and that it is the same as a transformation from an HGS frame with a
    # spherical representation

    obstime = "2011-01-01"
    hgscoord_cart = SkyCoord(x=1*u.km, y=0.*u.km, z=0.*u.km,
                             frame=HeliographicStonyhurst(obstime=obstime),
                             representation_type='cartesian')
    hpc_frame = Helioprojective(observer='earth', obstime=obstime)
    hgscoord_sph = hgscoord_cart.copy()
    hgscoord_sph.representation_type = 'spherical'
    hpccoord_cart = hgscoord_cart.transform_to(hpc_frame)
    hpccoord_sph = hgscoord_sph.transform_to(hpc_frame)
    assert_quantity_allclose(hpccoord_cart.Tx, hpccoord_sph.Tx)
    assert_quantity_allclose(hpccoord_cart.Ty, hpccoord_sph.Ty)
    assert_quantity_allclose(hpccoord_cart.distance, hpccoord_sph.distance)


def test_hgs_cartesian_rep_to_hcc():
    # This test checks transformation HGS->HCC when the coordinate is in a Cartesian
    # representation and that it is the same as a transformation from an HGS frame with a
    # spherical representation

    obstime = "2011-01-01"
    hgscoord_cart = SkyCoord(x=1*u.km, y=0.*u.km, z=0.*u.km,
                             frame=HeliographicStonyhurst(obstime=obstime),
                             representation_type='cartesian')
    hcc_frame = Heliocentric(observer='earth', obstime=obstime)
    hgscoord_sph = hgscoord_cart.copy()
    hgscoord_sph.representation_type = 'spherical'
    hcccoord_cart = hgscoord_cart.transform_to(hcc_frame)
    hcccoord_sph = hgscoord_sph.transform_to(hcc_frame)
    assert_quantity_allclose(hcccoord_cart.x, hcccoord_sph.x)
    assert_quantity_allclose(hcccoord_cart.y, hcccoord_sph.y)
    assert_quantity_allclose(hcccoord_cart.z, hcccoord_sph.z)


def test_hgs_cartesian_rep_to_hgc():
    # This test checks transformation HGS->HCC when the coordinate is in a Cartesian
    # representation and that it is the same as a transformation from an HGS frame with a
    # spherical representation

    obstime = "2011-01-01"
    hgscoord_cart = SkyCoord(x=1*u.km, y=0.*u.km, z=0.*u.km,
                             frame=HeliographicStonyhurst(obstime=obstime),
                             representation_type='cartesian')
    hgscoord_sph = hgscoord_cart.copy()
    hgscoord_sph.representation_type = 'spherical'
    # HGC
    hgcframe = HeliographicCarrington(observer='earth', obstime=obstime)
    hgccoord_cart = hgscoord_cart.transform_to(hgcframe)
    hgccoord_sph = hgscoord_sph.transform_to(hgcframe)
    assert_quantity_allclose(hgccoord_cart.lat, hgccoord_sph.lat)
    assert_quantity_allclose(hgccoord_cart.lon, hgccoord_sph.lon)
    assert_quantity_allclose(hgccoord_cart.radius, hgccoord_sph.radius)


def test_hcc_to_hpc_different_observer():
    # This test checks transformation HCC->HPC in the case where the HCC and HPC frames are
    # defined by different observers.

    rsun = 1*u.m
    D0 = 1*u.km
    L0 = 1*u.deg
    observer_1 = HeliographicStonyhurst(lat=0*u.deg, lon=0*u.deg, radius=D0)
    observer_2 = HeliographicStonyhurst(lat=0*u.deg, lon=L0, radius=D0)
    hcc_frame = Heliocentric(observer=observer_1)
    hpc_frame = Helioprojective(observer=observer_2)
    hcccoord = SkyCoord(x=rsun, y=rsun, z=rsun, frame=hcc_frame)
    hpccoord_out = hcccoord.transform_to(hpc_frame)
    hpccoord_expected = hcccoord.transform_to(HeliographicStonyhurst).transform_to(hpc_frame)
    assert_quantity_allclose(hpccoord_out.Tx, hpccoord_expected.Tx)
    assert_quantity_allclose(hpccoord_out.Ty, hpccoord_expected.Ty)
    assert_quantity_allclose(hpccoord_out.distance, hpccoord_expected.distance)


def test_hpc_to_hcc_different_observer():
    # This test checks transformation HPC->HCC in the case where the HCC and HPC frames are
    # defined by different observers.

    rsun = 1*u.m
    D0 = 1*u.km
    L0 = 1*u.deg
    observer_1 = HeliographicStonyhurst(lat=0*u.deg, lon=0*u.deg, radius=D0)
    observer_2 = HeliographicStonyhurst(lat=0*u.deg, lon=L0, radius=D0)
    hcc_frame = Heliocentric(observer=observer_1)
    hpc_frame = Helioprojective(observer=observer_2, rsun=rsun)
    hpccoord = SkyCoord(Tx=0*u.arcsec, Ty=0*u.arcsec, frame=hpc_frame)
    hcccoord_out = hpccoord.transform_to(hcc_frame)
    hcccoord_expected = hpccoord.transform_to(HeliographicStonyhurst).transform_to(hcc_frame)
    assert_quantity_allclose(hcccoord_out.x, hcccoord_expected.x)
    assert_quantity_allclose(hcccoord_out.y, hcccoord_expected.y)
    assert_quantity_allclose(hcccoord_out.z, hcccoord_expected.z)


def test_hcc_to_hpc_same_observer():
    # This test checks transformation HCC->HPC in the case of same observer

    rsun = 1*u.m
    D0 = 1*u.km
    observer = HeliographicStonyhurst(lat=0*u.deg, lon=0*u.deg, radius=D0)
    hcc_frame = Heliocentric(observer=observer)
    hpc_frame = Helioprojective(observer=observer, rsun=rsun)
    hcccoord = SkyCoord(x=rsun, y=rsun, z=rsun, frame=hcc_frame)
    hpccoord_out = hcccoord.transform_to(hpc_frame)
    hpccoord_expected = hcccoord.transform_to(HeliographicStonyhurst).transform_to(hpc_frame)
    assert_quantity_allclose(hpccoord_out.Tx, hpccoord_expected.Tx)
    assert_quantity_allclose(hpccoord_out.Ty, hpccoord_expected.Ty)
    assert_quantity_allclose(hpccoord_out.distance, hpccoord_expected.distance)


def test_hpc_to_hcc_same_observer():
    # This test checks transformation HPC->HCC in the case of same observer

    rsun = 1*u.m
    D0 = 1 * u.km
    observer = HeliographicStonyhurst(lat=0 * u.deg, lon=0 * u.deg, radius=D0)
    hcc_frame = Heliocentric(observer=observer)
    hpc_frame = Helioprojective(observer=observer, rsun=rsun)
    hpccoord = SkyCoord(Tx=0 * u.arcsec, Ty=0 * u.arcsec, frame=hpc_frame)
    hcccoord_out = hpccoord.transform_to(hcc_frame)
    hcccoord_expected = hpccoord.transform_to(HeliographicStonyhurst).transform_to(hcc_frame)
    assert_quantity_allclose(hcccoord_out.x, hcccoord_expected.x)
    assert_quantity_allclose(hcccoord_out.y, hcccoord_expected.y)
    assert_quantity_allclose(hcccoord_out.z, hcccoord_expected.z)


def test_hpc_hcc_different_observer_radius():
    # Tests HPC->HCC with a change in observer at different distances from the Sun
    observer1 = HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU)
    hpc = Helioprojective(0*u.arcsec, 0*u.arcsec, 0.5*u.AU, observer=observer1)

    observer2 = HeliographicStonyhurst(90*u.deg, 0*u.deg, 0.75*u.AU)
    hcc = hpc.transform_to(Heliocentric(observer=observer2))

    assert_quantity_allclose(hcc.x, -0.5*u.AU)
    assert_quantity_allclose(hcc.y, 0*u.AU, atol=1e-10*u.AU)
    assert_quantity_allclose(hcc.z, 0*u.AU, atol=1e-10*u.AU)


def test_hgs_hgs():
    # Test HGS loopback transformation
    obstime = Time('2001-01-01')
    old = SkyCoord(90*u.deg, 10*u.deg, 1*u.AU, frame=HeliographicStonyhurst(obstime=obstime))
    new = old.transform_to(HeliographicStonyhurst(obstime=obstime + 1*u.day))

    assert_quantity_allclose(new.lon, old.lon - 1*u.deg, atol=0.1*u.deg)  # due to Earth motion
    assert_quantity_allclose(new.lat, old.lat, atol=1e-3*u.deg)
    assert_quantity_allclose(new.radius, old.radius, atol=1e-5*u.AU)


def test_hgc_hgc():
    # Test HGC loopback transformation
    obstime = Time('2001-01-01')
    old = SkyCoord(90*u.deg, 10*u.deg, 1*u.AU, frame=HeliographicCarrington(observer='earth',
                                                                            obstime=obstime))
    new = old.transform_to(HeliographicCarrington(observer='earth', obstime=obstime + 1*u.day))

    assert_quantity_allclose(new.lon, 75.815607 * u.deg, atol=1e-7*u.deg)  # solar rotation
    # These are not equal to the old values, because the coordinates stay fixed
    # in inertial space, whilst the frame (fixed to the center of the Sun)
    # moves slightly.
    assert_quantity_allclose(new.lat, 9.999963 * u.deg, atol=1e-7*u.deg)
    assert_quantity_allclose(new.radius, 1.000009 * u.AU, atol=1e-7*u.AU)


def test_hgc_hgc_different_observers():
    obstime = Time('2001-01-01')
    hgc_earth = HeliographicCarrington(observer='earth', obstime=obstime)
    hgc_mars = HeliographicCarrington(observer='mars', obstime=obstime)
    hgc_sun = HeliographicCarrington(observer='sun', obstime=obstime)

    sc = SkyCoord(10*u.deg, 20*u.deg, 1*u.AU, frame=HeliographicStonyhurst(obstime=obstime))
    sc_hgc_earth = sc.transform_to(hgc_earth)

    sc_hgc_mars = sc_hgc_earth.transform_to(hgc_mars)

    sc_hgc_sun = sc_hgc_mars.transform_to(hgc_sun)

    ltt_earth = hgc_earth.observer.radius / speed_of_light
    assert_quantity_allclose(sc_hgc_earth.lon - sc_hgc_sun.lon, ltt_earth * sidereal_rotation_rate)

    ltt_mars = hgc_mars.observer.radius / speed_of_light
    assert_quantity_allclose(sc_hgc_mars.lon - sc_hgc_sun.lon, ltt_mars * sidereal_rotation_rate)


def test_hgc_self_observer():
    # Test specifying observer='self' for HGC
    obstime = Time('2001-01-01')
    hgc = HeliographicCarrington(10*u.deg, 20*u.deg, 3*u.AU, observer='self', obstime=obstime)

    # Transform to HGS (i.e., observer='self' in the source frame)
    hgs = hgc.transform_to(HeliographicStonyhurst(obstime=obstime))

    # Manually calculate the post-transformation longitude
    lon = sun.L0(obstime,
                 light_travel_time_correction=False,
                 nearest_point=False,
                 aberration_correction=False)
    lon += (hgc.radius - _RSUN) / speed_of_light * sidereal_rotation_rate

    assert_quantity_allclose(Longitude(hgs.lon + lon), hgc.lon)
    assert_quantity_allclose(hgs.lat, hgc.lat)
    assert_quantity_allclose(hgs.radius, hgc.radius)

    # Transform back to HGC (i.e., observer='self' in the destination frame)
    hgc_loop = hgs.transform_to(hgc.replicate_without_data())

    assert_quantity_allclose(hgc_loop.lon, hgc.lon)
    assert_quantity_allclose(hgc_loop.lat, hgc.lat)
    assert_quantity_allclose(hgc_loop.radius, hgc.radius)


def test_hgc_loopback_self_observer():
    # Test the HGC loopback where only one end has observer='self'
    obstime = Time('2001-01-01')
    coord = HeliographicCarrington(10*u.deg, 20*u.deg, 3*u.AU, observer='self', obstime=obstime)

    new_observer = HeliographicStonyhurst(40*u.deg, 50*u.deg, 6*u.AU)
    new_frame = HeliographicCarrington(observer=new_observer, obstime=obstime)

    new_coord = coord.transform_to(new_frame)

    # Manually calculate the longitude shift due to the difference in Sun-observer distance
    lon = (6*u.AU - 3*u.AU) / speed_of_light * sidereal_rotation_rate

    assert_quantity_allclose(new_coord.lon, coord.lon + lon)
    assert_quantity_allclose(new_coord.lat, coord.lat)
    assert_quantity_allclose(new_coord.radius, coord.radius)


def test_hcc_hcc():
    # Test same observer and changing obstime
    observer = HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU, obstime='2001-02-01')
    from_hcc = Heliocentric(0.2*u.AU, 0.3*u.AU, 0.4*u.AU, observer=observer, obstime='2001-01-01')
    to_hcc = from_hcc.transform_to(Heliocentric(observer=observer, obstime='2001-03-31'))

    # Since the observer is the same, the coordinates should be nearly the same but not exactly
    # equal due to motion of the origin (the Sun)
    assert np.all(from_hcc.cartesian.xyz != to_hcc.cartesian.xyz)
    assert_quantity_allclose(from_hcc.cartesian.xyz, to_hcc.cartesian.xyz, rtol=2e-3)

    # Test changing observer and same obstime
    observer1 = HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU, obstime='2001-01-01')
    observer2 = HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU, obstime='2001-03-31')
    from_hcc = Heliocentric(0.2*u.AU, 0.3*u.AU, 0.4*u.AU, observer=observer1, obstime='2001-02-01')
    to_hcc = from_hcc.transform_to(Heliocentric(observer=observer2, obstime='2001-02-01'))

    # This change in observer is approximately a 90-degree rotation about the Y axis
    assert_quantity_allclose(to_hcc.x, -from_hcc.z, rtol=2e-3)
    assert_quantity_allclose(to_hcc.y, from_hcc.y, rtol=2e-3)
    assert_quantity_allclose(to_hcc.z, from_hcc.x, rtol=2e-3)


def test_hcc_hgs_observer_mismatch():
    # Test whether the transformation gives the same answer regardless of what obstime the observer
    # coordinate is represented in
    observer1 = HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU, obstime='2001-01-01')
    observer2 = observer1.transform_to(HeliographicStonyhurst(obstime='2001-03-31'))

    hcc1 = Heliocentric(0.2*u.AU, 0.3*u.AU, 0.4*u.AU, observer=observer1, obstime=observer1.obstime)
    hgs1 = hcc1.transform_to(HeliographicStonyhurst(obstime=hcc1.obstime))

    hcc2 = Heliocentric(0.2*u.AU, 0.3*u.AU, 0.4*u.AU, observer=observer2, obstime=observer1.obstime)
    hgs2 = hcc2.transform_to(HeliographicStonyhurst(obstime=hcc2.obstime))

    assert_quantity_allclose(hgs1.lon, hgs2.lon)
    assert_quantity_allclose(hgs1.lat, hgs2.lat)
    assert_quantity_allclose(hgs1.radius, hgs2.radius)


def test_hgs_hcc_observer_mismatch():
    # Test whether the transformation gives the same answer regardless of what obstime the observer
    # coordinate is represented in
    observer1 = HeliographicStonyhurst(0*u.deg, 0*u.deg, 1*u.AU, obstime='2001-01-01')
    observer2 = observer1.transform_to(HeliographicStonyhurst(obstime='2001-03-31'))

    hgs = HeliographicStonyhurst(20*u.deg, 40*u.deg, 0.5*u.AU, obstime=observer1.obstime)
    hcc1 = hgs.transform_to(Heliocentric(observer=observer1, obstime=hgs.obstime))
    hcc2 = hgs.transform_to(Heliocentric(observer=observer2, obstime=hgs.obstime))

    assert_quantity_allclose(hcc1.cartesian.xyz, hcc2.cartesian.xyz)


def test_hgs_hcrs_sunspice():
    # Compare our HGS->HCRS transformation against SunSPICE by transforming beyond it
    # "HEQ" is another name for HEEQ, which is equivalent to Heliographic Stonyhurst
    # "HAE" is equivalent to Astropy's Heliocentric Mean Ecliptic, and defaults to J2000.0
    #
    # IDL> coord = [1.d, 0.d, 10.d]
    # IDL> convert_sunspice_lonlat, '2019-06-01', coord, 'HEQ', 'HAE', /au, /degrees
    # IDL> print, coord
    #        1.0000000      -108.65371       10.642778

    old = SkyCoord(0*u.deg, 10*u.deg, 1*u.AU, frame=HeliographicStonyhurst(obstime='2019-06-01'))
    new = old.transform_to(HeliocentricMeanEcliptic)

    assert_quantity_allclose(new.lon, Longitude(-108.65371*u.deg), atol=0.1*u.arcsec, rtol=0)
    assert_quantity_allclose(new.lat, 10.642778*u.deg, atol=0.1*u.arcsec, rtol=0)
    assert_quantity_allclose(new.distance, old.radius)

    # Transform to HAE precessed to the mean ecliptic of date instead of J2000.0
    # IDL> coord = [1.d, 0.d, 10.d]
    # IDL> convert_sunspice_lonlat, '2019-06-01', coord, 'HEQ', 'HAE', /precess, /au, /degrees
    # IDL> print, coord
    #        1.0000000      -108.38240       10.640314

    new = old.transform_to(HeliocentricMeanEcliptic(equinox='2019-06-01'))

    assert_quantity_allclose(new.lon, Longitude(-108.38240*u.deg), atol=0.1*u.arcsec, rtol=0)
    assert_quantity_allclose(new.lat, 10.640314*u.deg, atol=0.1*u.arcsec, rtol=0)
    assert_quantity_allclose(new.distance, old.radius)


def test_hgs_hgc_sunspice():
    # Compare our HGS->HGC transformation against SunSPICE
    # "HEQ" is another name for HEEQ, which is equivalent to Heliographic Stonyhurst
    # "Carrington" does not include light travel time to the observer, which our HGC includes
    #
    # IDL> coord = [1.d, 0.d, 10.d]
    # IDL> convert_sunspice_lonlat, '2019-06-01', coord, 'HEQ', 'Carrington', /au, /degrees
    # IDL> print, coord
    #        1.0000000       16.688242       10.000000

    old = SkyCoord(0*u.deg, 10*u.deg, 1*u.AU, frame=HeliographicStonyhurst(obstime='2019-06-01'))
    new = old.transform_to(HeliographicCarrington(observer='earth'))

    # Calculate the difference in longitude due to light travel time from the Sun to the Earth
    delta_lon = sidereal_rotation_rate * (sun.earth_distance(old.obstime) - _RSUN) / speed_of_light

    assert_quantity_allclose(new.lon, 16.688242*u.deg + delta_lon, atol=1e-2*u.arcsec, rtol=0)
    assert_quantity_allclose(new.lat, old.lat)
    assert_quantity_allclose(new.radius, old.radius)


def test_hgs_hcc_sunspice():
    # Compare our HGS->HCC transformation against SunSPICE
    # "HEQ" is another name for HEEQ, which is equivalent to Heliographic Stonyhurst
    # "HGRTN" is equivalent to our Heliocentric, but with the axes permuted
    # SunSPICE, like us, assumes an Earth observer if not explicitly specified
    #
    # IDL> coord = [7d5, 8d5, 9d5]
    # IDL> convert_sunspice_coord, '2019-06-01', coord, 'HEQ', 'HGRTN'
    # Assuming Earth observation
    # IDL> print, coord
    #        688539.32       800000.00       908797.89

    old = SkyCoord(CartesianRepresentation([7e5, 8e5, 9e5]*u.km),
                   frame=HeliographicStonyhurst(obstime='2019-06-01'))
    new = old.transform_to(Heliocentric(observer='earth'))

    assert_quantity_allclose(new.x, 800000.00*u.km, atol=1e-2*u.km)
    assert_quantity_allclose(new.y, 908797.89*u.km, atol=1e-2*u.km)
    assert_quantity_allclose(new.z, 688539.32*u.km, atol=1e-2*u.km)


def test_hpc_hgs_implicit_hcc():
    # An HPC->HGS transformation should give the same answer whether the transformation step
    #   through HCC is implicit or explicit
    start = SkyCoord(0*u.arcsec, 0*u.arcsec, 0.5*u.AU,
                     frame=Helioprojective(obstime='2019-06-01', observer='earth'))
    frame = HeliographicStonyhurst(obstime='2019-12-01')

    implicit = start.transform_to(frame)
    explicit1 = start.transform_to(Heliocentric(obstime=start.obstime, observer='earth')).\
        transform_to(frame)
    explicit2 = start.transform_to(Heliocentric(obstime=frame.obstime, observer='earth')).\
        transform_to(frame)

    assert_quantity_allclose(implicit.separation_3d(explicit1), 0*u.AU, atol=1e-10*u.AU)
    assert_quantity_allclose(implicit.separation_3d(explicit2), 0*u.AU, atol=1e-10*u.AU)


def test_velocity_hcrs_hgs():
    # Obtain the position/velocity of Earth in ICRS
    obstime = Time(['2019-01-01', '2019-04-01', '2019-07-01', '2019-10-01'])
    pos, vel = get_body_barycentric_posvel('earth', obstime)
    loc = pos.with_differentials(vel.represent_as(CartesianDifferential))
    earth = SkyCoord(loc, frame='icrs', obstime=obstime)

    # The velocity of Earth in HGS Y should be very close to zero because the XZ plane tracks
    # the Earth.
    new = earth.heliographic_stonyhurst
    assert_quantity_allclose(new.velocity.d_y, 0*u.km/u.s, atol=1e-5*u.km/u.s)

    # Test the loopback to ICRS
    newer = new.icrs
    assert_quantity_allclose(newer.velocity.d_xyz, vel.xyz)


def test_velocity_hgs_hgc():
    # Construct a simple HGS coordinate with zero velocity
    obstime = Time(['2019-01-01', '2019-04-01', '2019-07-01', '2019-10-01'])
    pos = CartesianRepresentation(1, 0, 0)*u.AU
    vel = CartesianDifferential(0, 0, 0)*u.km/u.s
    loc = (pos.with_differentials(vel))._apply('repeat', obstime.size)
    coord = SkyCoord(HeliographicStonyhurst(loc, obstime=obstime))

    # The induced velocity in HGC should be entirely longitudinal, and approximately equal to one
    # full rotation every mean synodic period (27.2753 days)
    hgc_frame = HeliographicCarrington(observer='earth', obstime=obstime)
    new = coord.transform_to(hgc_frame)
    new_vel = new.data.differentials['s'].represent_as(SphericalDifferential, new.data)
    assert_quantity_allclose(new_vel.d_lon, -360*u.deg / (27.27253*u.day), rtol=1e-2)
    assert_quantity_allclose(new_vel.d_lat, 0*u.deg/u.s)
    assert_quantity_allclose(new_vel.d_distance, 0*u.km/u.s, atol=1e-7*u.km/u.s)


def test_velocity_hgs_hci():
    # HGS and HCI share the same origin and Z axis, so the induced velocity is entirely angular
    obstime = Time(['2021-01-01', '2021-04-01', '2021-07-01', '2021-10-01'])
    venus_hgs = get_body_heliographic_stonyhurst('venus', obstime, include_velocity=True)
    venus_hci = venus_hgs.transform_to(HeliocentricInertial(obstime=obstime))

    # The induced velocity is the longitude component of Earth's velocity, ~360 deg/yr
    induced_dlon = get_earth(obstime, include_velocity=True).heliocentricinertial.d_lon
    assert_quantity_allclose(induced_dlon, 360*u.deg/u.yr, rtol=0.05)

    # The HCI velocity should be the same as the HGS velocity except for the induced velocity
    assert_quantity_allclose(venus_hci.d_distance, venus_hgs.d_radius, rtol=1e-5)
    assert_quantity_allclose(venus_hci.d_lon, venus_hgs.d_lon + induced_dlon, rtol=1e-6)
    assert_quantity_allclose(venus_hci.d_lat, venus_hgs.d_lat)


def test_velocity_hcrs_hci():
    # HCRS and HCI are both inertial frames with the same origin, so there is no induced velocity.
    # There is an induced angular velocity for HCRS->HGS, which should be canceled out by the
    # induced angular velocity for HGS->HCI.

    # Define an HCRS coordinate with a purely radial velocity
    sc_hcrs = SkyCoord(ra=[0, 90, 180, 270]*u.deg, pm_ra_cosdec=[0, 0, 0, 0]*u.deg/u.d,
                       dec=[10, 20, 30, 40]*u.deg, pm_dec=[0, 0, 0, 0]*u.deg/u.d,
                       distance=[5, 6, 7, 8]*u.AU, radial_velocity=[1, 2, 3, 4]*u.km/u.s,
                       frame='hcrs', obstime='2021-01-01')

    # Transform to HCI, and get the velocity vector in spherical coordinates
    sc_hci = sc_hcrs.heliocentricinertial

    # The HCI velocity should have the same amplitude, and should be purely radial
    assert_quantity_allclose(sc_hci.d_distance, sc_hcrs.velocity.norm(), rtol=1e-6)
    assert_quantity_allclose(sc_hci.d_lon, 0*u.arcsec/u.s, atol=1e-9*u.arcsec/u.s)
    assert_quantity_allclose(sc_hci.d_lat, 0*u.arcsec/u.s, atol=1e-9*u.arcsec/u.s)


def test_hme_hee_sunspice():
    # Compare our HME->HEE transformation against SunSPICE
    # "HAE" is equivalent to Astropy's Heliocentric Mean Ecliptic, and defaults to J2000.0
    #
    # IDL> coord = [1.d, 0.d, 10.d]
    # IDL> convert_sunspice_lonlat, '2019-06-01', coord, 'HAE', 'HEE', /au, /degrees
    # IDL> print, coord
    #        1.0000000       110.01610       10.000300

    old = SkyCoord(0*u.deg, 10*u.deg, 1*u.AU, frame=HeliocentricMeanEcliptic(obstime='2019-06-01'))
    new = old.transform_to(HeliocentricEarthEcliptic)

    assert_quantity_allclose(new.lon, Longitude(110.01610*u.deg), atol=0.01*u.arcsec, rtol=0)
    assert_quantity_allclose(new.lat, 10.000300*u.deg, atol=0.01*u.arcsec, rtol=0)
    assert_quantity_allclose(new.distance, old.distance)

    # Transform from HAE precessed to the mean ecliptic of date instead of J2000.0
    # IDL> coord = [1.d, 0.d, 10.d]
    # IDL> convert_sunspice_lonlat, '2019-06-01', coord, 'HAE', 'HEE', /au, /degrees, /precess
    # IDL> print, coord
    #        1.0000000       109.74535       10.000070

    old = SkyCoord(0*u.deg, 10*u.deg, 1*u.AU, frame=HeliocentricMeanEcliptic(obstime='2019-06-01',
                                                                             equinox='2019-06-01'))
    new = old.transform_to(HeliocentricEarthEcliptic)

    assert_quantity_allclose(new.lon, Longitude(109.74535*u.deg), atol=0.05*u.arcsec, rtol=0)
    assert_quantity_allclose(new.lat, 10.000070*u.deg, atol=0.01*u.arcsec, rtol=0)
    assert_quantity_allclose(new.distance, old.distance)


def test_hee_hee():
    # Test HEE loopback transformation
    obstime = Time('2001-01-01')
    old = SkyCoord(90*u.deg, 10*u.deg, 1*u.AU, frame=HeliocentricEarthEcliptic(obstime=obstime))

    new = old.transform_to(HeliocentricEarthEcliptic)

    assert_quantity_allclose(new.lon, old.lon)
    assert_quantity_allclose(new.lat, old.lat)
    assert_quantity_allclose(new.distance, old.distance)

    new = old.transform_to(HeliocentricEarthEcliptic(obstime=obstime + 1*u.day))

    assert_quantity_allclose(new.lon, old.lon - 1*u.deg, atol=0.1*u.deg)  # due to Earth motion
    assert_quantity_allclose(new.lat, old.lat, atol=0.5*u.arcsec)
    assert_quantity_allclose(new.distance, old.distance, rtol=1e-5)


def test_hee_gse_sunspice():
    # Compare our HEE->GSE transformation against SunSPICE
    #
    # IDL> coord = [0.7d, -20.d, 10.d]
    # IDL> convert_sunspice_coord, '2019-06-01', coord, 'HEE', 'GSE', /au, /degrees
    # IDL> print, coord
    #       0.45215884       32.777377       15.594639

    old = SkyCoord(-20*u.deg, 10*u.deg, 0.7*u.AU,
                   frame=HeliocentricEarthEcliptic(obstime='2019-06-01'))
    new = old.geocentricsolarecliptic

    assert_quantity_allclose(new.lon, 32.777377*u.deg, atol=0.01*u.arcsec, rtol=0)
    assert_quantity_allclose(new.lat, 15.594639*u.deg, atol=0.01*u.arcsec, rtol=0)
    assert_quantity_allclose(new.distance, 0.45215884*u.AU)


def test_gse_gse():
    # Test GSE loopback transformation
    old = SkyCoord(90*u.deg, 10*u.deg, 0.7*u.AU,
                   frame=GeocentricSolarEcliptic(obstime='2001-01-01'))
    new = old.transform_to(GeocentricSolarEcliptic)

    assert_quantity_allclose(new.lon, old.lon)
    assert_quantity_allclose(new.lat, old.lat)
    assert_quantity_allclose(new.distance, old.distance)


def test_hgs_hci_sunspice():
    # Compare our HGS->HCI transformation against SunSPICE
    # "HEQ" is another name for HEEQ, which is equivalent to Heliographic Stonyhurst
    #
    # IDL> coord = [1.d, 120.d, 10.d]
    # IDL> convert_sunspice_lonlat, '2019-06-01', coord, 'HEQ', 'HCI', /au, /degrees
    # IDL> print, coord
    #        1.0000000      -65.736793       10.000000

    old = SkyCoord(120*u.deg, 10*u.deg, 1*u.AU, frame=HeliographicStonyhurst(obstime='2019-06-01'))
    new = old.transform_to(HeliocentricInertial)

    assert_quantity_allclose(new.lon, -65.736793*u.deg, atol=0.5*u.arcsec, rtol=0)
    assert_quantity_allclose(new.lat, old.lat)
    assert_quantity_allclose(new.distance, old.radius)


def test_hci_hci():
    # Test HCI loopback transformation
    obstime = Time('2001-01-01')
    old = SkyCoord(90*u.deg, 10*u.deg, 0.7*u.AU, frame=HeliocentricInertial(obstime=obstime))
    new = old.transform_to(HeliocentricInertial)

    assert_quantity_allclose(new.lon, old.lon)
    assert_quantity_allclose(new.lat, old.lat)
    assert_quantity_allclose(new.distance, old.distance)

    new = old.transform_to(HeliocentricInertial(obstime=obstime + 1*u.day))

    assert_quantity_allclose(new.lon, old.lon, atol=0.1*u.deg)  # due to Earth motion
    assert_quantity_allclose(new.lat, old.lat, atol=1e-3*u.deg)
    assert_quantity_allclose(new.distance, old.distance, atol=1e-5*u.AU)


def test_hme_gei_sunspice():
    # Compare our HME->GEI transformation against SunSPICE
    # "HAE" is equivalent to Astropy's Heliocentric Mean Ecliptic, and defaults to J2000.0
    #
    # IDL> coord = [1.d, 120.d, 10.d]
    # IDL> convert_sunspice_lonlat, '2019-06-01', coord, 'HAE', 'GEI', /au, /degrees
    # IDL> print, coord
    #        1.8197210       95.230617       28.830109

    old = SkyCoord(120*u.deg, 10*u.deg, 1*u.AU,
                   frame=HeliocentricMeanEcliptic(obstime='2019-06-01'))
    new = old.transform_to(GeocentricEarthEquatorial)

    assert_quantity_allclose(new.lon, Longitude(95.230617*u.deg), atol=0.01*u.arcsec, rtol=0)
    assert_quantity_allclose(new.lat, 28.830109*u.deg, atol=0.05*u.arcsec, rtol=0)
    assert_quantity_allclose(new.distance, 1.8197210*u.AU)

    # Transform from HAE precessed to the mean ecliptic of date instead of J2000.0
    # IDL> coord = [1.d, 120.d, 10.d]
    # IDL> convert_sunspice_lonlat, '2019-06-01', coord, 'HAE', 'GEI', /au, /degrees, /precess
    # IDL> print, coord
    #        1.8217103       95.079030       28.827750

    old = SkyCoord(120*u.deg, 10*u.deg, 1*u.AU,
                   frame=HeliocentricMeanEcliptic(obstime='2019-06-01', equinox='2019-06-01'))
    new = old.transform_to(GeocentricEarthEquatorial(equinox=_J2000))

    assert_quantity_allclose(new.lon, Longitude(95.079030*u.deg), atol=0.05*u.arcsec, rtol=0)
    assert_quantity_allclose(new.lat, 28.827750*u.deg, atol=0.05*u.arcsec, rtol=0)
    assert_quantity_allclose(new.distance, 1.8217103*u.AU)


def test_gei_gei():
    # Test GEI loopback transformation using the 2017 revision to Franz & Harper 2002
    t = Time('1996-08-28 16:46:00', scale='tt')
    gei_j2000 = CartesianRepresentation([-5.7840451, -4.1082375, 1.9146822] * (6378.14*u.km))
    gei_d = CartesianRepresentation([-5.7864918, -4.1039136, 1.9165612] * (6378.14*u.km))

    old = SkyCoord(gei_j2000, frame=GeocentricEarthEquatorial(obstime=t))
    new = old.transform_to(GeocentricEarthEquatorial(equinox=t, obstime=t)).cartesian

    assert_quantity_allclose(new.xyz, gei_d.xyz)


def test_no_observer():
    # Tests transformations to and from observer-based frames with no observer defined
    frames_in = [Heliocentric(0*u.km, 0*u.km, 0*u.km, observer=None),
                 Heliocentric(0*u.km, 0*u.km, 0*u.km, observer=None, obstime='2001-01-01'),
                 Helioprojective(0*u.deg, 0*u.deg, observer=None),
                 Helioprojective(0*u.deg, 0*u.deg, observer=None, obstime='2001-01-01')]
    frames_out = frames_in + [
        HeliographicStonyhurst(0*u.deg, 0*u.deg, obstime=None),
        HeliographicStonyhurst(0*u.deg, 0*u.deg, obstime='2001-01-01'),
        Heliocentric(0*u.km, 0*u.km, 0*u.km, observer=None, obstime='2012-12-12'),
        Heliocentric(0*u.km, 0*u.km, 0*u.km, observer="earth", obstime=None),
        Heliocentric(0*u.km, 0*u.km, 0*u.km, observer="earth", obstime='2001-01-01'),
        Helioprojective(0*u.deg, 0*u.deg, observer=None, obstime='2012-12-12'),
        Helioprojective(0*u.deg, 0*u.deg, observer="earth", obstime=None),
        Helioprojective(0*u.deg, 0*u.deg, observer="earth", obstime='2001-01-01')]

    # Self-transformations should succeed
    for f in frames_in:
        f.transform_to(f.replicate_without_data())

    # All other transformations should error
    for i, f1 in enumerate(frames_in):
        for f2 in frames_out[i + 1:]:
            with pytest.raises(ConvertError):
                f1.transform_to(f2)
            with pytest.raises(ConvertError):
                f2.transform_to(f1)


def test_array_obstime():
    # Validate that you can transform from an array of obstimes to no obstimes,
    # or different obstimes.
    a = SkyCoord([10]*2, [10]*2, unit=u.deg,
                 observer="earth",
                 obstime=["2019-01-01", "2019-01-02"],
                 frame="heliographic_carrington")

    t = a.transform_to(Helioprojective)
    assert isinstance(t.frame, Helioprojective)

    t2 = a.transform_to(Helioprojective(obstime=["2019-01-03", "2019-01-04"]))
    assert isinstance(t2.frame, Helioprojective)


_frames_wo_observer = [HeliographicStonyhurst, HeliocentricInertial,
                       HeliocentricEarthEcliptic, GeocentricSolarEcliptic,
                       GeocentricEarthEquatorial]


@pytest.mark.parametrize("frame_class", _frames_wo_observer)
def test_convert_error_with_no_obstime(frame_class):
    # For most transformations, we do not allow `obstime` to be `None`
    frame = frame_class(CartesianRepresentation(0, 0, 0)*u.km, obstime=None)

    with pytest.raises(ConvertError, match=r".*obstime.*"):
        ICRS(0*u.deg, 0*u.deg, 0*u.AU).transform_to(frame)

    with pytest.raises(ConvertError, match=r".*obstime.*"):
        frame.transform_to(ICRS())


# Convenience function to check whether a transformation succeeds if the target `obstime` is `None`
def assert_no_obstime_on_target_end(start_class, end_class):
    start_obstime = Time("2001-01-01")

    if hasattr(start_class, 'observer'):
        coord = start_class(CartesianRepresentation(0, 0, 0)*u.km,
                            obstime=start_obstime, observer="earth")
    else:
        coord = start_class(CartesianRepresentation(0, 0, 0)*u.km, obstime=start_obstime)

    result = coord.transform_to(end_class(obstime=None))
    assert result.obstime == start_obstime


# We currently allow the target `obstime` to be `None` for the transformation subgraph
# below `HeliographicStonyhurst`, but this may change in the future
_frameset1 = [HeliographicStonyhurst, HeliocentricInertial]
_frameset2 = [HeliographicCarrington, Heliocentric, Helioprojective]


@pytest.mark.parametrize("start_class", _frameset1 + _frameset2)
@pytest.mark.parametrize("end_class", _frameset1)
def test_no_obstime_on_target_end_hgs_subgraph(start_class, end_class):
    assert_no_obstime_on_target_end(start_class, end_class)


# We currently allow the target `obstime` to be `None` for the transformation subgraph
# below `HeliocentricEarthEcliptic`, but this may change in the future
_frameset3 = [HeliocentricEarthEcliptic, GeocentricSolarEcliptic]


@pytest.mark.parametrize("start_class", _frameset3)
@pytest.mark.parametrize("end_class", _frameset3)
def test_no_obstime_on_target_end_hee_subgraph(start_class, end_class):
    assert_no_obstime_on_target_end(start_class, end_class)


def test_transform_with_sun_center():
    sun_center = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU,
                          frame=HeliographicStonyhurst(obstime="2001-01-01"))

    with transform_with_sun_center():
        result1 = sun_center.transform_to(HeliographicStonyhurst(obstime="2001-02-01"))

    # The coordinate should stay pointing at Sun center
    assert_quantity_allclose(result1.lon, sun_center.lon)
    assert_quantity_allclose(result1.lat, sun_center.lat)
    assert_quantity_allclose(result1.radius, sun_center.radius)

    other = SkyCoord(10*u.deg, 20*u.deg, 1*u.AU,
                     frame=HeliographicStonyhurst(obstime="2001-01-01"))

    with transform_with_sun_center():
        result2 = other.transform_to(HeliographicCarrington(observer='earth', obstime="2001-02-01"))

    # The coordinate should stay at the same latitude and the same distance from Sun center
    assert_quantity_allclose(result2.lat, other.lat)
    assert_quantity_allclose(result2.radius, other.radius)


def test_transform_with_sun_center_reset():
    # This test sequence ensures that the context manager resets propoerly

    sun_center = SkyCoord(0*u.deg, 0*u.deg, 0*u.AU,
                          frame=HeliographicStonyhurst(obstime="2001-01-01"))
    end_frame = HeliocentricInertial(obstime="2001-02-01")

    # Without the context manager, the coordinate should not point at Sun center
    result1 = sun_center.transform_to(end_frame)
    assert result1.lon != sun_center.lon
    assert result1.lat != sun_center.lat
    assert result1.distance != sun_center.radius

    # Using the context manager, the coordinate should point at Sun center
    with transform_with_sun_center():
        result2 = sun_center.transform_to(end_frame)
    assert_quantity_allclose(result2.lon, sun_center.lon)
    assert_quantity_allclose(result2.lat, sun_center.lat)
    assert_quantity_allclose(result2.distance, sun_center.radius)

    # Exiting a nested context manager should not affect the outer context manager
    with transform_with_sun_center():
        with transform_with_sun_center():
            pass
        result2a = sun_center.transform_to(end_frame)
    assert_quantity_allclose(result2a.lon, result2.lon)
    assert_quantity_allclose(result2a.lat, result2.lat)
    assert_quantity_allclose(result2a.distance, result2.distance)

    # After the context manager, the coordinate should have the same result as the first transform
    result3 = sun_center.transform_to(end_frame)
    assert_quantity_allclose(result3.lon, result1.lon)
    assert_quantity_allclose(result3.lat, result1.lat)
    assert_quantity_allclose(result3.distance, result1.distance)
