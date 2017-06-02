import numpy as np

import astropy.units as u
from astropy.tests.helper import quantity_allclose

from sunpy.coordinates import Helioprojective, HeliographicStonyhurst
from astropy.coordinates import SkyCoord


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

    sc_in = SkyCoord(0*u.arcsec, 0*u.arcsec, rsun=rsun, observer=observer_in, frame='helioprojective')
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
