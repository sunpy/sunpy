import numpy as np

import astropy.units as u
from astropy.tests.helper import quantity_allclose, assert_quantity_allclose
from astropy.coordinates import SkyCoord, get_body_barycentric
from astropy.time import Time

from sunpy.coordinates import (Helioprojective, HelioprojectiveRadial,
                               HeliographicStonyhurst, HeliographicCarrington,
                               get_sun_L0)
from sunpy.time import parse_time


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


def test_hgs_hgc_roundtrip():
    obstime = "2011-01-01"

    hgsin = HeliographicStonyhurst(lat=0*u.deg, lon=0*u.deg, obstime=obstime)
    hgcout = hgsin.transform_to(HeliographicCarrington)

    assert_quantity_allclose(hgsin.lat, hgcout.lat)
    assert_quantity_allclose(hgcout.lon, get_sun_L0(obstime))

    hgsout = hgcout.transform_to(HeliographicStonyhurst)

    assert_quantity_allclose(hgsout.lat, hgsin.lat)
    assert_quantity_allclose(hgsout.lon, hgsin.lon)


def test_hpc_hpr_roundtrip():
    hpc = Helioprojective(0*u.arcsec, -100*u.arcsec)

    hpr = hpc.transform_to(HelioprojectiveRadial)

    hpc2 = hpr.transform_to(Helioprojective)

    assert_quantity_allclose(hpc.Tx, hpc2.Tx, atol=1e-13*u.arcsec)
    assert_quantity_allclose(hpc.Ty, hpc2.Ty, atol=1e-13*u.arcsec)


def test_hpc_hpr_roundtrip_3d():
    hpc = Helioprojective(0*u.arcsec, -100*u.arcsec, obstime="2017-12-25")
    hpc = hpc.calculate_distance()

    hpr = hpc.transform_to(HelioprojectiveRadial)

    hpc2 = hpr.transform_to(Helioprojective)

    assert_quantity_allclose(hpc.Tx, hpc2.Tx, atol=1e-13*u.arcsec)
    assert_quantity_allclose(hpc.Ty, hpc2.Ty, atol=1e-13*u.arcsec)
    assert_quantity_allclose(hpc.distance, hpc2.distance, atol=1e-13*u.km)


def test_hpr_calculate_distance():
    hpc = Helioprojective(0*u.arcsec, -100*u.arcsec, obstime="2017-12-25")
    hpcd = hpc.calculate_distance()

    hpr = hpc.transform_to(HelioprojectiveRadial(obstime="2017-12-25"))
    hpr = hpr.calculate_distance()

    assert_quantity_allclose(hpcd.distance, hpr.distance, rtol=1e-4)
