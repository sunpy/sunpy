"""
Generic LightCurve Tests
"""
from __future__ import absolute_import

#
# @TODO:
#   ndarray + indices?
#   time deltas instead of datetimes?

#pylint: disable=C0103,R0904,W0201,W0232,E1101,E1103
import sunpy
import pytest
import datetime

def test_input_empty():
    with pytest.raises(NotImplementedError):
        sunpy.lightcurve.LightCurve()

def test_input_dict_datetimes():
    """Tests LightCurve creation from a dictionary and list of datetimes"""
    base = datetime.datetime.today()
    dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    sunpy.lightcurve.LightCurve({"param": range(24 * 60)}, index=dates)
    
        
        