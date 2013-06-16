"""
EVE tests
"""
from __future__ import absolute_import

#pylint: disable=C0103,R0904,W0201,W0232,E1103
import sunpy
from sunpy.data.test import (EVE_AVERAGES_CSV)
import matplotlib


def test_eve():
    eve = sunpy.lightcurve.EVELightCurve.create('2013/04/15')
    assert isinstance(eve, sunpy.lightcurve.EVELightCurve)
   
def test_txt():
    """Check support for parsing EVE TXT files """
    eve = sunpy.lightcurve.EVELightCurve.create(
    "http://lasp.colorado.edu/eve/data_access/quicklook/quicklook_data/L0CS/LATEST_EVE_L0CS_DIODES_1m.txt") 
    assert isinstance(eve, sunpy.lightcurve.EVELightCurve)        

def test_csv_parsing():
    """Check support for parsing EVE CSV files"""
    csv = sunpy.lightcurve.EVELightCurve.create(EVE_AVERAGES_CSV)
    assert isinstance(csv, sunpy.lightcurve.sources.eve.EVELightCurve)
