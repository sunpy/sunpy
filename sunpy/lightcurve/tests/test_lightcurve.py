"""
Generic LightCurve Tests
"""
from __future__ import absolute_import

#
# @TODO:
#   time deltas instead of datetimes?

#pylint: disable=C0103,R0904,W0201,W0232,E1101,E1103
import numpy as np
import pytest
import datetime
import sunpy
import sunpy.lightcurve
from sunpy.data.test import (EVE_AVERAGES_CSV)
import pandas

# Generate input test data
base = datetime.datetime.today()
dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]


@pytest.mark.parametrize(("data", "index"), [
    (range(24 * 60), dates),
    (np.arange(24 * 60), dates),
    ({"param": range(24 * 60)}, dates),
    ({"instr1": range(24*60), "instr2": range(24*60), "instr3": range(24*60)}, dates),
    ([{'instr1': x, 'instr2': x+1, 'instr3': x+2} for x in range(24*60)], dates),
    (pandas.Series(range(24*60)), dates),

])
def test_input(data, index):
    """Tests different types of expected input"""
    lc = sunpy.lightcurve.LightCurve.create(data, index=index)
    assert isinstance(lc, sunpy.lightcurve.LightCurve)
    assert isinstance(lc.data, pandas.DataFrame)
    assert len(lc.data.index) == 24*60
    assert lc.data.index[0] == base
    assert lc.data.index[-1] == base - datetime.timedelta(minutes=24*60-1)


@pytest.mark.parametrize(("bad_input"), [
    (None),
    (EVE_AVERAGES_CSV)
])
def test_unimplemented(bad_input):
    """Tests input that has not been implemented for the generic LC class"""
    with pytest.raises((TypeError, NotImplementedError)):
        sunpy.lightcurve.LightCurve.create(bad_input)

