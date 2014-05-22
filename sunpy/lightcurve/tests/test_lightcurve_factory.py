# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 18:29:53 2014

@author: stuart
"""
import os

import pytest

import sunpy.data.test
import sunpy.time
from sunpy.lightcurve.lightcurve_factory import LightCurve
from sunpy.lightcurve.sources import *
from sunpy.util.odict import OrderedDict
import glob
from astropy.io import fits

base_path = os.path.join(sunpy.data.test.rootdir, 'lightcurve/')

#def test_generic_dhp():
#    hdus = fits.open(os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'))
#    data = hdus[0].data
#    header = OrderedDict(hdus[0].header)
#    header.pop('INSTRUME', None)
#    
#    lc = LightCurve(data,header)
#    assert isinstance(lc, LYRALightCurve)
#
#def test_lyra_dhp():
#    hdus = fits.open(os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'))
#    data = hdus[0].data
#    header = OrderedDict(hdus[0].header)
#    
#    lc = LightCurve((data,header))
#    assert isinstance(lc, LYRALightCurve)

def test_goes_single():
    lc = LightCurve(os.path.join(base_path,'g15_xrs_2s_20120101_20120101_shrunk.csv'), source='goes')

    assert isinstance(lc, GOESLightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2012-01-01 00:17:11.547000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-01 23:54:09.617000')

def test_eve_single_as_genericlc():
    lc = LightCurve(os.path.join(base_path,'EVE_He_II_304_averages.csv'))

    assert isinstance(lc, GenericLightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2012-06-13T00:00:00Z')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-06-13T01:39:00Z')
    

def test_goes_webapi():
    lc = LightCurve(os.path.join(base_path,'xrs_2s_webapi_shrunk.csv'), source='goes')

    assert isinstance(lc, GOESLightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2012-06-01 00:17:03')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-06-02 23:50:07')

def test_lyra_single():
    lc = LightCurve(os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'))

    assert isinstance(lc, LYRALightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2012-01-01 00:00:00.118000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-01 23:59:27.164000')

def test_eve_single():
    lc = LightCurve(os.path.join(base_path, 'LATEST_EVE_L0CS_DIODES_1m.txt'), source='eve')

    assert isinstance(lc, EVELightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2012-06-19 00:00:00')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-06-19 01:40:00')

def test_norh_single():
    lc = LightCurve(os.path.join(base_path, 'norh_tcz_120101.fits'))

    assert isinstance(lc, NoRHLightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2011-12-31 22:44:50.641000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-01 06:30:01.641000')
    
def test_norh_single_34():
    lc = LightCurve(os.path.join(base_path, 'norh_tca_120101.fits'))

    assert isinstance(lc, NoRHLightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2011-12-31 22:44:50.591000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-01 06:30:01.591000')

def test_multi_source():
    lc = LightCurve(os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                    os.path.join(base_path,'g15_xrs_2s_20120101_20120101_shrunk.csv'),
                                      source=['lyra', 'goes'])
    assert isinstance(lc[0], LYRALightCurve)
    assert isinstance(lc[1], GOESLightCurve)

def test_multi():
    lc = LightCurve(os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                    os.path.join(base_path,'g15_xrs_2s_20120101_20120101_shrunk.csv'),
                                      source=[None, 'goes'])
    assert isinstance(lc[0], LYRALightCurve)
    assert isinstance(lc[1], GOESLightCurve)

def give_source_for_fname(fname):
    if 'eve' in fname:
    	return 'eve'
    if 'lyra' in fname:
    	return 'lyra'
    if 'norh' in fname: 
    	return 'norh'
    if 'xrs' in fname:
    	return 'goes'
    return None


@pytest.mark.parametrize(("inp", "source"),
[('*.fits', None),
 #('*.csv', [give_source_for_fname(x.lower()) for x in glob.glob(base_path + '*.csv')])]) eve_csv file not recognised
 ('*.csv', None),
 ('', None ]) # test for directories
def test_glob_dir(inp,source):
    oeve = EVELightCurve
    ogoes = GOESLightCurve
    onorh = NoRHLightCurve
    olyra = LYRALightCurve
    map_ = {}
    map_['eve'] = oeve
    map_['geos'] = ogoes
    map_['norh'] = onorh
    map_['lyra'] = olyra
    lc = LightCurve(base_path + inp,source=source)
    #for i,ilc in enumerate(lc):                     lc is not iterable, concatenate flag required 
    #    assert isinstance(ilc,map_[source[i]])


def test_one_from_many():
    lc = LightCurve([os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                     os.path.join(base_path,'lyra_20120102-000000_lev2_std.fits')])

    assert isinstance(lc, LYRALightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2012-01-01T00:00:00.118000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-02 23:59:52.188000')

def test_many_from_many():
    lc = sunpy.lightcurve.LightCurve([os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                                      os.path.join(base_path,'lyra_20120102-000000_lev2_std.fits')],
                                     [os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                                      os.path.join(base_path,'lyra_20120102-000000_lev2_std.fits')])

    assert isinstance(lc[0], LYRALightCurve)
    assert isinstance(lc[1], LYRALightCurve)
    assert lc[0].time_range().start() == sunpy.time.parse_time('2012-01-01T00:00:00.118000')
    assert lc[0].time_range().end() == sunpy.time.parse_time('2012-01-02 23:59:52.188000')
    assert lc[1].time_range().start() == sunpy.time.parse_time('2012-01-01T00:00:00.118000')
    assert lc[1].time_range().end() == sunpy.time.parse_time('2012-01-02 23:59:52.188000')


#* Create a LightCurve from a file and then truncate to a time range.
def truncate_file():
    lc = sunpy.lightcurve.LightCurve(os.path.join(base_path,'lyra_20120101-000000_lev2_std.fits'),
                                      timerange=["2012/1/1T00:00:00", "2012/1/1T12:00:00"])

    assert isinstance(lc, LYRALightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2012-01-01 00:00:00.118000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-01 11:59:59.964000')

@pytest.mark.online
def test_source_date():
    lc = LightCurve(timerange='2012/1/2', source='lyra')

    assert isinstance(lc, LYRALightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2012-01-02 00:00:00.110000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-02 23:59:58.988000')

@pytest.mark.online
def test_get_norh_34():
    lc = LightCurve(timerange='2012/1/1', wavelength='34', source='norh')

    assert isinstance(lc, NoRHLightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2011-12-31 22:44:50.641000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-01 05:14:51.641000')

@pytest.mark.online
def test_get_norh():
    lc = LightCurve(timerange='2012/1/1', source='norh')

    assert isinstance(lc, NoRHLightCurve)
    assert lc.time_range().start() == sunpy.time.parse_time('2011-12-31 22:44:50.591000')
    assert lc.time_range().end() == sunpy.time.parse_time('2012-01-01 05:14:51.591000')
