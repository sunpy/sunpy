# -*- coding: utf-8 -*-

from astropy.wcs import WCS

from ..frames import *
from ..wcs_utils import solar_wcs_frame_mapping


def test_hpr():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['HRLN', 'HRLT']

    result = solar_wcs_frame_mapping(wcs)

    assert isinstance(result, HelioprojectiveRadial)


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
