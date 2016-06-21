import numpy as np

import astropy.units as u
from astropy.tests.helper import quantity_allclose

from sunpy.coordinates import Helioprojective


def test_hpc_hpc():
    # Use some unphysical values for solar parameters for testing
    rsun = 1*u.m
    D0 = 1*u.km
    L0 = 1*u.deg

    hpc_in = Helioprojective(0*u.arcsec, 0*u.arcsec, rsun=rsun, D0=D0)
    hpc_out = Helioprojective(L0=L0, D0=D0, rsun=rsun)

    hpc_new = hpc_in.transform_to(hpc_out)

    assert hpc_new.L0 == hpc_out.L0
    assert hpc_new.B0 == hpc_out.B0
    assert hpc_new.D0 == hpc_out.D0

    # Calculate the distance subtended by an angle of L0 from the centre of the
    # Sun.
    dd = -1 * rsun * np.tan(L0)
    # Calculate the angle corresponding to that distance as seen by the new
    # observer.
    theta = np.arctan2(dd, (D0 - rsun))

    assert quantity_allclose(theta, hpc_new.Tx, rtol=1e-3)


def test_hpc_hpc_null():

    hpc_in = Helioprojective(0*u.arcsec, 0*u.arcsec)
    hpc_out = Helioprojective()

    hpc_new = hpc_in.transform_to(hpc_out)

    assert hpc_new is not hpc_in
    assert quantity_allclose(hpc_new.Tx, hpc_in.Tx)
    assert quantity_allclose(hpc_new.Ty, hpc_in.Ty)
    assert quantity_allclose(hpc_new.D0, hpc_in.D0)
    assert quantity_allclose(hpc_new.B0, hpc_in.B0)
    assert quantity_allclose(hpc_new.L0, hpc_in.L0)
