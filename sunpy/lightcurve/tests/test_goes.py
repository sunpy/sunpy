
"""
GOES LightCurve Tests
"""
from __future__ import absolute_import

import pytest
import sunpy.lightcurve
import sunpy.time

class TestGOESLightCurve():
    
    @pytest.mark.online
    def test_goes_range(self):
       pytest.skip('The GOES API currently does not work')
       lc1 = sunpy.lightcurve.GOESLightCurve.create('2011/06/01', '2011/06/02')
       
       assert isinstance(lc1, sunpy.lightcurve.GOESLightCurve)

    @pytest.mark.online
    def test_goes_timerange(self):
        pytest.skip('The GOES API currently does not work')
        timerange = sunpy.time.TimeRange('2011/06/01', '2011/06/02')
        lc1 = sunpy.lightcurve.GOESLightCurve.create(timerange)
       
        assert isinstance(lc1, sunpy.lightcurve.GOESLightCurve)
    
    def compare(self,lc1,lc2):
        try:
            (lc1.data == lc2.data)
        except:
            raise Exception

    @pytest.mark.online
    def test_filename(self):
        pytest.skip('The GOES API currently does not work')
        lc1 = sunpy.lightcurve.GOESLightCurve.create('2011/06/01', '2011/06/02')
        lc2 = sunpy.lightcurve.GOESLightCurve.create('2011/06/03', '2011/06/04')
        #If the dataframes are non-idential it raises an error, if they are
        #identical it returns True
        with pytest.raises((Exception)):
            self.compare(lc1, lc2)
        
        #This snippet is better but python 2.7 only, maybe sometime it will be 
        #useful
#        with self.assertRaises(Exception) as context:
#            (lc1.data == lc2.data).all()
#        
#        self.assertEqual(context.exception.message,
#                         'Can only compare identically-labeled DataFrame objects')
