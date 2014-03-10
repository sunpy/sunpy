# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 18:29:53 2014

@author: stuart
"""

import os
import sunpy.time
from sunpy.lightcurve.lightcurve_factory import LightCurve
from sunpy.lightcurve.sources import *

base_path = '/home/stuart/sunpy/data/'

def test_goes_single():
    lc = LightCurve(os.path.join(base_path,'xrs_2s.csv'), source='goes')
    assert isinstance(lc, GOESLightCurve)

def test_lyra_single():
    lc = LightCurve(os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'))

    assert isinstance(lc, LYRALightCurve)

def test_eve_single():
    lc = LightCurve(os.path.join(base_path, 'LATEST_EVE_L0CS_DIODES_1m.txt'), source='eve')

    assert isinstance(lc, EVELightCurve)

def test_norh_single():
    lc = LightCurve(os.path.join(base_path, 'norh_tca120101_test.fits'))

    assert isinstance(lc, NoRHLightCurve)

def test_multi_source():
    lc = LightCurve(os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                    os.path.join(base_path,'g15_xrs_2s_20120101_20120102.csv'),
                                      source=['lyra', 'goes'])
    assert isinstance(lc[0], LYRALightCurve)
    assert isinstance(lc[1], GOESLightCurve)

def test_multi():
    lc = LightCurve(os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                    os.path.join(base_path,'g15_xrs_2s_20120101_20120102.csv'),
                                      source=[None, 'goes'])
    assert isinstance(lc[0], LYRALightCurve)
    assert isinstance(lc[1], GOESLightCurve)

def test_one_from_many():
    lc = LightCurve([os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                     os.path.join(base_path,'lyra_20120102-000000_lev2_std.fits')])

    assert isinstance(lc, LYRALightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2012-01-01T00:00:00.118000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-02 23:59:59.988000')

def test_many_from_many():
    lc = sunpy.lightcurve.LightCurve([os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                                      os.path.join(base_path,'lyra_20120102-000000_lev2_std.fits')],
                                     [os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                                      os.path.join(base_path,'lyra_20120102-000000_lev2_std.fits')])

    assert isinstance(lc[0], LYRALightCurve)
    assert isinstance(lc[1], LYRALightCurve)
    assert lc[0].time_range().start() == sunpy.time.parse_time('2012-01-01T00:00:00.118000')
    assert lc[0].time_range().end() == sunpy.time.parse_time('2012-01-02 23:59:59.988000')
    assert lc[1].time_range().start() == sunpy.time.parse_time('2012-01-01T00:00:00.118000')
    assert lc[1].time_range().end() == sunpy.time.parse_time('2012-01-02 23:59:59.988000')


#* Create a LightCurve from a file and then truncate to a time range.
def truncate_file():
    lc = sunpy.lightcurve.LightCurve(os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                                      timerange=["2012/1/1T00:00:00", "2012/1/1T12:00:00"])

    assert isinstance(lc, LYRALightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2012-01-01 00:00:00.118000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-01 11:59:59.964000')




def test_source_date():
    lc = LightCurve(timerange='2012/1/2', source='goes')

    assert isinstance(lc, GOESLightCurve)

def test_get_norh_34():
    lc = LightCurve(timerange='2012/1/1', wavelength='34', source='norh')

    assert isinstance(lc, NoRHLightCurve)