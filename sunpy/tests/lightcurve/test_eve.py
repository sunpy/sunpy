"""
EVE tests
"""
from __future__ import absolute_import

#pylint: disable=C0103,R0904,W0201,W0232,E1103
import sunpy
from sunpy.data.test import (EVE_AVERAGES_CSV)

class TestEve:
    """Tests the EVE class"""
    def setup_class(self):
        pass

    def teardown_class(self):
        pass
        
    def test_csv_parsing(self):
        """Check support for parsing EVE CSV files"""
        sunpy.lightcurve.EVELightCurve.create(EVE_AVERAGES_CSV)
