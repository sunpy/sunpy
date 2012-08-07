"""
Generic LightCurve Tests
"""
from __future__ import absolute_import

#
# @TODO:
#   ndarray + indices?
#   time deltas instead of datetimes?

#pylint: disable=C0103,R0904,W0201,W0232,E1101,E1103
import numpy as np
import pytest
import datetime
import sunpy
import pandas
from sunpy.data.test import (EVE_AVERAGES_CSV)

def test_input_empty():
    """Tests empty input"""
    with pytest.raises(NotImplementedError):
        sunpy.lightcurve.LightCurve()

def test_input_dict_datetimes():
    """Tests LightCurve creation from a dictionary and list of datetimes"""
    base = datetime.datetime.today()
    dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    sunpy.lightcurve.LightCurve({"param": range(24 * 60)}, index=dates)
    
def test_input_ndarray_datetimes():
    """Tests LightCurve creation from a dictionary and list of datetimes"""
    base = datetime.datetime.today()
    dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    sunpy.lightcurve.LightCurve(np.arange(24 * 60), index=dates)
    
def test_input_list_datetimes():
    """Tests LightCurve creation from a dictionary and list of datetimes"""
    base = datetime.datetime.today()
    dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    sunpy.lightcurve.LightCurve(range(24 * 60), index=dates)
    
def test_input_file():
    """Tests filepath input"""
    with pytest.raises(NotImplementedError):
        sunpy.lightcurve.LightCurve(EVE_AVERAGES_CSV)
        
