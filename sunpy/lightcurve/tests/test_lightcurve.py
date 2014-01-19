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

# Generate input test data
base = datetime.datetime.today()
dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]

@pytest.mark.parametrize(("data", "index"), [
    (range(24 * 60), dates),
    (np.arange(24 * 60), dates),
    ({"param": range(24 * 60)}, dates)
])
def test_input(data, index):
    """Tests different types of expected input"""
    sunpy.lightcurve.LightCurve.create(data, index=index)

@pytest.mark.parametrize(("bad_input"), [
    (None),
    (EVE_AVERAGES_CSV)
])
def test_unimplemented(bad_input):
    """Tests input that has not been implemented for the generic LC class"""
    with pytest.raises((TypeError, NotImplementedError)):
        sunpy.lightcurve.LightCurve.create(bad_input)
