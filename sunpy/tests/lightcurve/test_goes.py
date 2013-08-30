"""
GOES LightCurve Tests
"""
from __future__ import absolute_import

import sunpy.lightcurve

class TestGOESLightCurve():
    def test_goes():
        goes = sunpy.lightcurve.EVELightCurve.create('2013/04/15')
        assert isinstance(eve, sunpy.lightcurve.GOESLightCurve)
    
    def test_filename():
        goes = sunpy.lightcurve.EVELightCurve.create('2013/04/15')
        assert isinstance(eve, sunpy.lightcurve.GOESLightCurve)
