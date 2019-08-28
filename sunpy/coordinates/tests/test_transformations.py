import numpy as np
import pytest

import astropy
import astropy.units as u
from astropy.tests.helper import quantity_allclose, assert_quantity_allclose
from astropy.coordinates import (SkyCoord, get_body_barycentric, Angle,
                                 ConvertError, Longitude, CartesianRepresentation,
                                 get_body_barycentric_posvel,
                                 CartesianDifferential, SphericalDifferential)
# Versions of Astropy that do not have HeliocentricMeanEcliptic have the same frame
# with the misleading name HeliocentricTrueEcliptic
try:
    from astropy.coordinates import HeliocentricMeanEcliptic
except ImportError:
    from astropy.coordinates import HeliocentricTrueEcliptic as HeliocentricMeanEcliptic

from astropy.time import Time

from sunpy.coordinates import (Helioprojective, HeliographicStonyhurst,
                               HeliographicCarrington, Heliocentric,
                               get_earth)
from sunpy.coordinates import sun
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
    hgs_out = hcc_in.transform_to(HeliographicStonyhurst)

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
    hgcout = hgsin.transform_to(HeliographicCarrington(obstime=obstime))

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
    hgscoord_sph = hgscoord_cart.copy()
    hgscoord_sph.representation_type = 'spherical'
    hpccoord_cart = hgscoord_cart.transform_to(Helioprojective(obstime=obstime))
    hpccoord_sph = hgscoord_sph.transform_to(Helioprojective(obstime=obstime))
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
    hgscoord_sph = hgscoord_cart.copy()
    hgscoord_sph.representation_type = 'spherical'
    hcccoord_cart = hgscoord_cart.transform_to(Heliocentric(obstime=obstime))
    hcccoord_sph = hgscoord_sph.transform_to(Heliocentric(obstime=obstime))
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
    hgccoord_cart = hgscoord_cart.transform_to(HeliographicCarrington(obstime=obstime))
    hgccoord_sph = hgscoord_sph.transform_to(HeliographicCarrington(obstime=obstime))
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
    old = SkyCoord(90*u.deg, 10*u.deg, 1*u.AU, frame=HeliographicCarrington(obstime=obstime))
    new = old.transform_to(HeliographicCarrington(obstime=obstime + 1*u.day))

    assert_quantity_allclose(new.lon, old.lon - 14.1844*u.deg, atol=1e-4*u.deg)  # solar rotation
    assert_quantity_allclose(new.lat, old.lat, atol=1e-4*u.deg)
    assert_quantity_allclose(new.radius, old.radius, atol=1e-5*u.AU)


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
    # "Carrington" is offset by 0.076 degrees in longitude from our Heliographic Carrington (HGC)
    #   because "Carrington" does not include light travel time to the observer, while our
    #   HGC includes the light travel time to Earth (see Seidelmann et al. 2007).
    #
    # IDL> coord = [1.d, 0.d, 10.d]
    # IDL> convert_sunspice_lonlat, '2019-06-01', coord, 'HEQ', 'Carrington', /au, /degrees
    # IDL> print, coord
    #        1.0000000       16.688242       10.000000

    old = SkyCoord(0*u.deg, 10*u.deg, 1*u.AU, frame=HeliographicStonyhurst(obstime='2019-06-01'))
    new = old.heliographic_carrington

    assert_quantity_allclose(new.lon, 16.688242*u.deg + 0.076*u.deg, atol=1e-2*u.arcsec, rtol=0)
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
    new = old.heliocentric

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


@pytest.mark.skipif(astropy.__version__ < '3.2.0', reason="Not supported by Astropy <3.2")
def test_velocity_hcrs_hgs():
    # Obtain the position/velocity of Earth in ICRS
    obstime = Time(['2019-01-01', '2019-04-01', '2019-07-01', '2019-10-01'])
    pos, vel = get_body_barycentric_posvel('earth', obstime)
    loc = pos.with_differentials(vel.represent_as(CartesianDifferential))
    earth = SkyCoord(loc, frame='icrs', obstime=obstime)

    # The velocity of Earth in HGS should be very close to zero.  The velocity in the HGS Y
    # direction is slightly further away from zero because there is true latitudinal motion.
    new = earth.heliographic_stonyhurst
    assert_quantity_allclose(new.velocity.d_x, 0*u.km/u.s, atol=1e-15*u.km/u.s)
    assert_quantity_allclose(new.velocity.d_y, 0*u.km/u.s, atol=1e-14*u.km/u.s)
    assert_quantity_allclose(new.velocity.d_x, 0*u.km/u.s, atol=1e-15*u.km/u.s)

    # Test the loopback to ICRS
    newer = new.icrs
    assert_quantity_allclose(newer.velocity.d_x, vel.x)
    assert_quantity_allclose(newer.velocity.d_y, vel.y)
    assert_quantity_allclose(newer.velocity.d_z, vel.z)


def test_velocity_hgs_hgc():
    # Construct a simple HGS coordinate with zero velocity
    obstime = Time(['2019-01-01', '2019-04-01', '2019-07-01', '2019-10-01'])
    pos = CartesianRepresentation(1, 0, 0)*u.AU
    vel = CartesianDifferential(0, 0, 0)*u.km/u.s
    loc = (pos.with_differentials(vel))._apply('repeat', obstime.size)
    coord = SkyCoord(HeliographicStonyhurst(loc, obstime=obstime))

    # The induced velocity in HGC should be entirely longitudinal, and approximately equal to one
    # full rotation every mean synodic period (27.2753 days)
    new = coord.heliographic_carrington
    new_vel = new.data.differentials['s'].represent_as(SphericalDifferential, new.data)
    assert_quantity_allclose(new_vel.d_lon, -360*u.deg / (27.27253*u.day), rtol=1e-2)
    assert_quantity_allclose(new_vel.d_lat, 0*u.deg/u.s)
    assert_quantity_allclose(new_vel.d_distance, 0*u.km/u.s, atol=1e-7*u.km/u.s)
