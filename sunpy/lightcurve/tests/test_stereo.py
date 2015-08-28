"""
STEREO LightCurve Tests

"""

# This module was developed with funding from 
# Google Summer of Code 2015
# author - Ankit Kumar  <ankitkmr.iitk@gmail.com>

from __future__ import absolute_import

import pytest
import os

import sunpy.lightcurve
from sunpy.time import TimeRange
import sunpy.data.test

filepath = sunpy.data.test.rootdir


class TestLETLightCurve(object):

    @pytest.mark.parametrize("data_file",[('Al_narrow_ahead.txt'),('Ar_summed_ahead.txt'),
    ('Fe_sectored_ahead_2015_001_level1_11.txt'), ('Al_summed_ahead.txt'), ('Ar_ahead_2006_318_level1_11.txt'), 
    ('Ar_summed_ahead_2015_1hr_level1_11.txt'), ('Fe_sectored_ahead_2015_01_10min_level1_11.txt'),
    ('Al_summed_ahead_2007_01_10min_level1_11.txt'), ('Ar_narrow_ahead.txt'), ('CNO_lo_sectored_ahead_2015_1hr_level1_11.txt')])
    @pytest.mark.online
    def test_header(self, data_file):
        """Test header parsing from file"""
        lc = sunpy.lightcurve.LETLightCurve
        assert (lc._parse_txt(os.path.join(filepath ,'let' ,data_file))[0])[1][:14] in ['Flux for Bin 0', 'Column 6: LET '] 


    @pytest.mark.parametrize("data_file",[('Al_narrow_ahead.txt'),('Ar_summed_ahead.txt'),
    ('Fe_sectored_ahead_2015_001_level1_11.txt'), ('Al_summed_ahead.txt'), ('Ar_ahead_2006_318_level1_11.txt'), 
    ('Ar_summed_ahead_2015_1hr_level1_11.txt'), ('Fe_sectored_ahead_2015_01_10min_level1_11.txt'),
    ('Al_summed_ahead_2007_01_10min_level1_11.txt'), ('Ar_narrow_ahead.txt'), ('CNO_lo_sectored_ahead_2015_1hr_level1_11.txt')])
    @pytest.mark.online
    def test_data(self, data_file):
        """Test for non empty data parsing from file"""
        lc = sunpy.lightcurve.LETLightCurve
        assert (lc._parse_txt(os.path.join(filepath ,'let' ,data_file))[1]).empty == False



class TestSITLightCurve(object):

    @pytest.mark.parametrize("data_file",[('SIT_Ahead_10min_4HE_2007_01.txt'), ('SIT_Ahead_10min_Fe_2007_01.txt'),
    ('SIT_Ahead_10min_H_2007_01.txt'), ('SIT_Ahead_10min_O_2007_01.txt')])
    @pytest.mark.online
    def test_header(self, data_file):
        """Test parsed header and data columns list for equal lengths """
        lc = sunpy.lightcurve.SITLightCurve._parse_txt(os.path.join(filepath ,'sit' ,data_file))
        assert len(lc[0]) == len(lc[1].columns) #length of header list equals number of columns


    @pytest.mark.parametrize("data_file",[('SIT_Ahead_10min_4HE_2007_01.txt'), ('SIT_Ahead_10min_Fe_2007_01.txt'),
    ('SIT_Ahead_10min_H_2007_01.txt'), ('SIT_Ahead_10min_O_2007_01.txt')])    
    @pytest.mark.online
    def test_data(self, data_file):
        """Test for non empty data parsing from file"""
        lc = sunpy.lightcurve.SITLightCurve
        assert (lc._parse_txt(os.path.join(filepath ,'sit' ,data_file))[1]).empty == False


class TestHETLightCurve(object):

    @pytest.mark.parametrize("data_file",[('AeH06Dec.12h.txt'), ('AeH06Dec.15m.txt'),
    ('AeH06Dec.1d.txt'), ('AeH06Dec.1h.txt'), ('AeH06Dec.1m.txt')])
    @pytest.mark.online
    def test_header(self, data_file):
        """Test parsed header and data columns list for equal lengths """
        lc = sunpy.lightcurve.HETLightCurve._parse_txt(os.path.join(filepath ,'het' ,data_file))
        assert len(lc[0]) == len(lc[1].columns) #length of header list equals number of columns


    @pytest.mark.parametrize("data_file",[('AeH06Dec.12h.txt'), ('AeH06Dec.15m.txt'),
    ('AeH06Dec.1d.txt'), ('AeH06Dec.1h.txt'), ('AeH06Dec.1m.txt')])
    @pytest.mark.online
    def test_data(self, data_file):
        """Test for non empty data parsing from file"""
        lc = sunpy.lightcurve.HETLightCurve
        assert (lc._parse_txt(os.path.join(filepath ,'het' ,data_file))[1]).empty == False


class TestPLASTICLightCurve(object):

    @pytest.mark.parametrize("data_file",[('STA_L2_PLA_1DMax_10min_20140101_001_V09.txt'),  
    ('STA_L2_PLA_1DMax_1min_20140101_001_V09.txt'),('STA_L2_PLA_1DMax_1hr_20141016_289_V09.txt')])
    @pytest.mark.online
    def test_header(self, data_file):
        """Test header parsing from file"""
        lc = sunpy.lightcurve.PLASTICLightCurve
        assert lc._parse_txt(os.path.join(filepath ,'plastic' ,data_file))[0][-2:] == [ 'Spcrft. Long. [HCI, degrees]', 'Spcrft. Lat. [HCI, degrees]' ]


    @pytest.mark.parametrize("data_file",[('STA_L2_PLA_1DMax_10min_20140101_001_V09.txt'),  
    ('STA_L2_PLA_1DMax_1min_20140101_001_V09.txt'),('STA_L2_PLA_1DMax_1hr_20141016_289_V09.txt')])
    @pytest.mark.online
    def test_data(self, data_file):
        """Test for non empty data parsing from file"""
        lc = sunpy.lightcurve.PLASTICLightCurve
        assert (lc._parse_txt(os.path.join(filepath , 'plastic' ,data_file))[1]).empty == False


    @pytest.mark.parametrize("data_file",[('STA_L2_PLA_1DMax_10min_20140101_001_V09.txt'),  
    ('STA_L2_PLA_1DMax_1min_20140101_001_V09.txt'),('STA_L2_PLA_1DMax_1hr_20141016_289_V09.txt')])
    @pytest.mark.online
    def test_header(self, data_file):
        """Test parsed header and data columns list for equal lengths """
        lc = sunpy.lightcurve.PLASTICLightCurve._parse_txt(os.path.join(filepath ,'plastic' ,data_file))
        assert len(lc[0]) == len(lc[1].columns) #length of header list equals number of columns


class TestSEPTLightCurve(object):

    @pytest.mark.parametrize("data_file",[('sept_ahead_ele_asun_2015_001_10min_l2_v03.dat.txt'), ('sept_ahead_ele_asun_2015_001_1h_l2_v03.dat.txt'),
    ('sept_ahead_ele_asun_2015_001_1d_l2_v03.dat.txt'), ('sept_ahead_ele_asun_2015_001_1min_l2_v03.dat.txt')])
    @pytest.mark.online
    def test_data(self, data_file):
        """Test for non empty data parsing from file"""
        lc = sunpy.lightcurve.SEPTLightCurve
        assert (lc._parse_txt(os.path.join(filepath ,'sept' ,data_file))[1]).empty == False


    @pytest.mark.parametrize("data_file",[('sept_ahead_ele_asun_2015_001_10min_l2_v03.dat.txt'), ('sept_ahead_ele_asun_2015_001_1h_l2_v03.dat.txt'),
    ('sept_ahead_ele_asun_2015_001_1d_l2_v03.dat.txt'), ('sept_ahead_ele_asun_2015_001_1min_l2_v03.dat.txt')])
    @pytest.mark.online
    def test_header(self, data_file):
        """Test parsed header and data columns list for equal lengths """
        lc = sunpy.lightcurve.SEPTLightCurve._parse_txt(os.path.join(filepath ,'sept' ,data_file))
        assert len(lc[0]) == len(lc[1].columns) #length of header list equals number of columns


