# -*- coding: utf-8 -*-

from astropy.wcs import WCS

from ..frames import *
from ..wcs_utils import solar_wcs_frame_mapping

def test_hpc():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['HPLN', 'HPLT']

    result = solar_wcs_frame_mapping(wcs)

    assert result is HelioProjective

def test_hgs():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['HGLN', 'HGLT']

    result = solar_wcs_frame_mapping(wcs)

    assert result is HelioGraphicStonyhurst

def test_hgc():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['CRLN', 'CRLT']

    result = solar_wcs_frame_mapping(wcs)

    assert result is HelioGraphicCarrington

def test_hcc():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['SOLX', 'SOLY']

    result = solar_wcs_frame_mapping(wcs)

    assert result is HelioCentric

def test_none():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['wibble', 'wobbl']

    result = solar_wcs_frame_mapping(wcs)

    assert result is None
