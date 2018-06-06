import numpy as np
import pytest

import astropy.units as u
from astropy.tests.helper import quantity_allclose, assert_quantity_allclose
from astropy.coordinates import SkyCoord, get_body_barycentric, HeliocentricTrueEcliptic, Angle
from astropy.time import Time

from sunpy.coordinates import (Helioprojective, HeliographicStonyhurst,
                               HeliographicCarrington, Heliocentric,
                               get_sun_L0, get_earth)
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

    assert hpc_new.observer == hpc_out.observer


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
    now = Time(parse_time('now'))
    earth_hcrs = SkyCoord(get_body_barycentric('earth', now), frame='icrs', obstime=now).hcrs

    # Convert from HCRS to HGS
    earth_hgs = earth_hcrs.transform_to(HeliographicStonyhurst)

    # The HGS longitude of the Earth should be zero within numerical error
    assert quantity_allclose(earth_hgs.lon, 0*u.deg, atol=1e-12*u.deg)

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
    # HeliocentricTrueEcliptic (HTE).  It will fail if there are errors in Astropy's
    # HCRS->ICRS or ICRS->HTE transformations.

    # Use published HGS coordinates in the Astronomical Almanac (2013), pages C6-C7
    obstime = Time('2013-01-28')
    earth_hgs = SkyCoord(0*u.deg, -5.73*u.deg, 0.9848139*u.AU, frame=HeliographicStonyhurst,
                         obstime=obstime)

    # Transform to HTE at observation-time equinox
    earth_hte = earth_hgs.transform_to(HeliocentricTrueEcliptic(equinox=obstime))

    # Validate against published values from the Astronomical Almanac (2013), page C6 per page E2
    # The dominant source of inaccuracy is the limited precision of the published B0 used above
    assert quantity_allclose(earth_hte.lon, Angle('308d13m30.51s') - 180*u.deg, atol=5*u.arcsec)
    assert quantity_allclose(earth_hte.lat, -Angle('-0.27s'), atol=10*u.arcsec)
    assert quantity_allclose(earth_hte.distance, 0.9848139*u.AU, atol=5e-7*u.AU)


def test_hgs_hgc_roundtrip():
    obstime = "2011-01-01"

    hgsin = HeliographicStonyhurst(lat=0*u.deg, lon=0*u.deg, obstime=obstime)
    hgcout = hgsin.transform_to(HeliographicCarrington(obstime=obstime))

    assert_quantity_allclose(hgsin.lat, hgcout.lat)
    assert_quantity_allclose(hgcout.lon, get_sun_L0(obstime))

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
                             representation='cartesian')
    hgscoord_sph = hgscoord_cart.copy()
    hgscoord_sph.representation = 'spherical'
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
                             representation='cartesian')
    hgscoord_sph = hgscoord_cart.copy()
    hgscoord_sph.representation = 'spherical'
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
                             representation='cartesian')
    hgscoord_sph = hgscoord_cart.copy()
    hgscoord_sph.representation = 'spherical'
    # HGC
    hgccoord_cart = hgscoord_cart.transform_to(HeliographicCarrington(obstime=obstime))
    hgccoord_sph = hgscoord_sph.transform_to(HeliographicCarrington(obstime=obstime))
    assert_quantity_allclose(hgccoord_cart.lat, hgccoord_sph.lat)
    assert_quantity_allclose(hgccoord_cart.lon, hgccoord_sph.lon)
    assert_quantity_allclose(hgccoord_cart.radius, hgccoord_sph.radius)
    