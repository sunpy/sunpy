from __future__ import absolute_import

import sunpy
from sunpy.roi import roi

def test_roi_instance():
    region=roi(times=['2012-06-20 05:00','2012-06-20 07:00'],description='dummy_roi')
    assert isinstance(region,sunpy.roi.roi)

    
