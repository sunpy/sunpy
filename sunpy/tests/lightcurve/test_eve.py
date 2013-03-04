"""
EVE tests
"""
from __future__ import absolute_import

#pylint: disable=C0103,R0904,W0201,W0232,E1103
import sunpy
import unittest
from sunpy.data.test import (EVE_AVERAGES_CSV)

class TestEve(unittest.TestCase):
    """Tests the EVE class"""
    
    def test_eve(self):
        eve = sunpy.lightcurve.EVELightCurve.create('2012/06/20')
        self.assertIsInstance(eve, sunpy.lightcurve.EVELightCurve)

    def test_txt(self):
        """Check support for parsing EVE TXT files """
        eve = sunpy.lightcurve.EVELightCurve.create(
        "http://lasp.colorado.edu/eve/data_access/quicklook/quicklook_data/L0CS/LATEST_EVE_L0CS_DIODES_1m.txt") 
        self.assertIsInstance(eve, sunpy.lightcurve.EVELightCurve)        
        
    def test_csv_parsing(self):
        """Check support for parsing EVE CSV files"""
        csv = sunpy.lightcurve.EVELightCurve.create(EVE_AVERAGES_CSV)
        self.assertIsInstance(csv, sunpy.lightcurve.sources.eve.EVELightCurve)