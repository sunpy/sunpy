"""
NOAA LightCurve Tests
"""
from __future__ import absolute_import

import pytest
import sunpy.lightcurve
import unittest

class TestNOAAIndicesLightCurve(unittest.TestCase):
    
    @pytest.mark.online
    def test_create(self):
       lc = sunpy.lightcurve.NOAAIndicesLightCurve.create()
       assert isinstance(lc, sunpy.lightcurve.NOAAIndicesLightCurve)
    
    @pytest.mark.online
    def test_isempty(self):
        lc = sunpy.lightcurve.NOAAIndicesLightCurve.create()
        self.assertFalse(lc.data.empty)
        
class TestNOAAPredictIndicesLightCurve(unittest.TestCase):

    @pytest.mark.online
    def test_create(self):
       lc = sunpy.lightcurve.NOAAPredictIndicesLightCurve.create()
       assert isinstance(lc, sunpy.lightcurve.NOAAPredictIndicesLightCurve)
    
    @pytest.mark.online
    def test_isempty(self):
        lc = sunpy.lightcurve.NOAAPredictIndicesLightCurve.create()
        self.assertFalse(lc.data.empty)
        