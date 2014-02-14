"""
GOES LightCurve Tests
"""
from __future__ import absolute_import

import pytest
import sunpy.lightcurve


class TestGOESLightCurve():
    
    @pytest.mark.online
    def test_goes_range(self):
       lc1 = sunpy.lightcurve.GOESLightCurve.create('2012/06/01','2012/06/02')
       
       assert isinstance(lc1, sunpy.lightcurve.GOESLightCurve)
    
    def compare(self,lc1,lc2):
        try:
            (lc1.data == lc2.data)
        except:
            raise Exception

    @pytest.mark.online
    def test_filename(self):
        lc1 = sunpy.lightcurve.GOESLightCurve.create('2012/06/01','2012/06/02')
        lc2 = sunpy.lightcurve.GOESLightCurve.create('2012/06/03','2012/06/04')
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
