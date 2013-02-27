"""
Lyra Tests
"""
from __future__ import absolute_import

#pylint: disable=C0103,R0904,W0201,W0232,E1103
import sunpy
import unittest
import matplotlib as mpl

class TestLyra(unittest.TestCase):
    """Tests the Lyra class"""
    
    def test_lyra(self):
        lyra = sunpy.lightcurve.LYRALightCurve.create(
        "http://proba2.oma.be/lyra/data/bsd/2011/08/10/lyra_20110810-000000_lev2_std.fits")
        self.assertIsInstance(lyra, sunpy.lightcurve.LYRALightCurve)
        
    def test_peek(self):
        lyra = sunpy.lightcurve.LYRALightCurve.create(
        "http://proba2.oma.be/lyra/data/bsd/2011/08/10/lyra_20110810-000000_lev2_std.fits")
        self.assertIsInstance(lyra.peek(),mpl.figure.Figure)