"""
This file is a stop-gap test file that compares this implementation to the old
implementation in sunpy.wcs.

The idea is that sunpy.wcs will be deprecated and removed, so this should not
be relied on for actual testing of the transformation framework.
"""

from ..frames import Helioprojective, Heliocentric, HeliographicStonyhurst
from sunpy import wcs
from sunpy.tests.helpers import assert_quantity_allclose

import astropy.units as u

import pytest


@pytest.mark.parametrize('Tx, Ty', [(0*u.arcsec, 0*u.arcsec),
                                    (10*u.arcsec, 0*u.arcsec),
                                    (-100*u.arcsec, -1000*u.arcsec),
                                    (40.0*u.arcsec, 32.0*u.arcsec),
                                    (1500*u.arcsec, 1500*u.arcsec)])
def test_hpc_hcc(Tx, Ty):
    hpc = Helioprojective(Tx, Ty)
    hcc = hpc.transform_to(Heliocentric)

    x, y, z = wcs.convert_hpc_hcc(Tx.value, Ty.value, angle_units='arcsec',
                                  dsun_meters=hpc.D0.to(u.m), z=True)

    assert_quantity_allclose(x*u.m, hcc.x)
    assert_quantity_allclose(y*u.m, hcc.y)
    assert_quantity_allclose(z*u.m, hcc.z)



@pytest.mark.parametrize('Tx, Ty', [(0*u.arcsec, 0*u.arcsec),
                                    (10*u.arcsec, 0*u.arcsec),
                                    (-100*u.arcsec, -1000*u.arcsec),
                                    (40.0*u.arcsec, 32.0*u.arcsec),
                                    (1500*u.arcsec, 1500*u.arcsec)])
def test_hpc_hgs(Tx, Ty):
    hpc = Helioprojective(Tx, Ty)
    hgs = hpc.transform_to(HeliographicStonyhurst)

    lon, lat = wcs.convert_hpc_hg(Tx.value, Ty.value, angle_units='arcsec',
                                  b0_deg=hpc.B0.to(u.deg).value, l0_deg=hpc.L0.to(u.deg).value,
                                  dsun_meters=hpc.D0.to(u.m))

    assert_quantity_allclose(lon*u.deg, hgs.lon)
    assert_quantity_allclose(lat*u.deg, hgs.lat)


@pytest.mark.parametrize('lon, lat', [(0*u.deg, 0*u.deg),
                                      (34.0*u.deg, 45.0*u.deg),
                                      (-80*u.deg, 70*u.deg)])
def test_hgs_hpc(lon, lat):
    hgs = HeliographicStonyhurst(lon, lat)
    hpc = hgs.transform_to(Helioprojective)

    Tx, Ty = wcs.convert_hg_hpc(lon.value, lat.value, angle_units='arcsec',
                                  b0_deg=hpc.B0.to(u.deg).value, l0_deg=hpc.L0.to(u.deg).value,
                                  dsun_meters=hpc.D0.to(u.m))

    assert_quantity_allclose(Tx*u.arcsec, hpc.Tx)
    assert_quantity_allclose(Ty*u.arcsec, hpc.Ty)


@pytest.mark.parametrize('lon, lat', [(0*u.deg, 0*u.deg),
                                      (34.0*u.deg, 45.0*u.deg),
                                      (-80*u.deg, 70*u.deg)])
def test_hgs_hcc(lon, lat):
    hgs = HeliographicStonyhurst(lon, lat)
    hcc = hgs.transform_to(Heliocentric)

    x, y, z = wcs.convert_hg_hcc(lon.value, lat.value,
                                 r=hgs.radius.to(u.m).value,
                                 z=True)

    assert_quantity_allclose(x*u.m, hcc.x)
    assert_quantity_allclose(y*u.m, hcc.y)
    assert_quantity_allclose(z*u.m, hcc.z)
