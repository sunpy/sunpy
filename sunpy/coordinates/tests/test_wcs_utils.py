# -*- coding: utf-8 -*-

import numpy as np

import sunpy.map
from astropy.wcs import WCS

from sunpy.coordinates.frames import Helioprojective, Heliocentric, HeliographicStonyhurst, HeliographicCarrington
from ..wcs_utils import solar_wcs_frame_mapping


def test_hpc():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['HPLN', 'HPLT']

    result = solar_wcs_frame_mapping(wcs)

    assert isinstance(result, Helioprojective)


def test_hgs():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['HGLN', 'HGLT']

    result = solar_wcs_frame_mapping(wcs)

    assert isinstance(result, HeliographicStonyhurst)


def test_hgc():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['CRLN', 'CRLT']

    result = solar_wcs_frame_mapping(wcs)

    assert isinstance(result, HeliographicCarrington)


def test_hcc():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['SOLX', 'SOLY']

    result = solar_wcs_frame_mapping(wcs)

    assert isinstance(result, Heliocentric)


def test_none():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['spam', 'eggs']

    result = solar_wcs_frame_mapping(wcs)

    assert result is None


def test_wcs_extras():
    """
    To enable proper creation of the coordinate systems, Map sticks three extra
    attributes on the WCS object:
    * heliographic_longitude
    * heliographic_latitude
    * dsun
    """
    data = np.ones([6, 6], dtype=np.float64)
    header = {'CRVAL1': 0,
              'CRVAL2': 0,
              'CRPIX1': 5,
              'CRPIX2': 5,
              'CDELT1': 10,
              'CDELT2': 10,
              'CUNIT1': 'arcsec',
              'CUNIT2': 'arcsec',
              'PC1_1': 0,
              'PC1_2': -1,
              'PC2_1': 1,
              'PC2_2': 0,
              'NAXIS1': 6,
              'NAXIS2': 6,
              'CTYPE1': 'HPLN-TAN',
              'CTYPE2': 'HPLT-TAN',
              'date-obs': '1970/01/01T00:00:00',
              'obsrvtry': 'Foo',
              'detector': 'bar',
              'wavelnth': 10,
              'waveunit': 'm',
              'hglt_obs': 0,
              'hgln_obs': 0,
              'dsun_obs': 10,
              'rsun_ref': 690000000}
    generic_map = sunpy.map.Map((data, header))

    wcs = generic_map.wcs

    assert wcs.heliographic_observer.lat.value == 0
    assert wcs.heliographic_observer.lon.value == 0
    assert wcs.heliographic_observer.radius.value == 10
    assert wcs.rsun.value == header['rsun_ref']

    result = solar_wcs_frame_mapping(wcs)

    assert isinstance(result, Helioprojective)
    assert result.observer.lat.value == 0
    assert result.observer.lon.value == 0
    assert result.observer.radius.value == 10
    assert result.rsun.value == header['rsun_ref']
