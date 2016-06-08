from sunpy.coordinates import Helioprojective

import astropy.units as u
from astropy.tests.helper import quantity_allclose

import numpy as np


def test_hpc_hpc():
    rsun = 1*u.m
    D0 = 1*u.km
    L0 = 1*u.deg

    hpc_in = Helioprojective(0*u.arcsec, 0*u.arcsec, rsun=rsun, D0=D0)
    hpc_out = Helioprojective(L0=1*u.deg, D0=D0, rsun=rsun)

    hpc_new = hpc_in.transform_to(hpc_out)

    assert hpc_new.L0 == hpc_out.L0
    assert hpc_new.B0 == hpc_out.B0
    assert hpc_new.D0 == hpc_out.D0

    # Calculate the distance subtended by an angle of L0 from the centre of the
    # Sun.
    dd = -1 * rsun * np.tan(L0)
    # Calculate the angle corresponding to that distance as seen by the new
    # observer.
    theta = np.arctan(dd / D0)

    assert quantity_allclose(theta, hpc_new.Tx, rtol=1e-3)
