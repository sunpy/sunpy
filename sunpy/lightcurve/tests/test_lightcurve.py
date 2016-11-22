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
dates = [base + datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
base_input = np.arange(24 * 60)

@pytest.mark.parametrize(("data", "index"), [
    (base_input, dates),
    (base_input, dates),
    ({"param": base_input}, dates),
    ({"instr1": base_input, "instr2": base_input, "instr3": base_input}, dates),
    ([{'instr1': x, 'instr2': x+1, 'instr3': x+2} for x in base_input], dates),
    (pandas.Series(base_input), dates),

])
def test_input(data, index):
    """Tests different types of expected input"""
    lc = sunpy.lightcurve.LightCurve.create(data, index=index)
    assert isinstance(lc, sunpy.lightcurve.LightCurve)
    assert isinstance(lc.data, pandas.DataFrame)
    assert len(lc.data.index) == 24 * 60
    assert lc.data.index[0] == base
    assert lc.data.index[-1] == base + datetime.timedelta(minutes=24 * 60 - 1)


@pytest.mark.parametrize(("bad_input"), [
    (None),
    (EVE_AVERAGES_CSV)
])
def test_unimplemented(bad_input):
    """Tests input that has not been implemented for the generic LC class"""
    with pytest.raises((TypeError, NotImplementedError)):
        sunpy.lightcurve.LightCurve.create(bad_input)


@pytest.mark.parametrize(("overlap_factor"), [0, 1, 5, 10, 100, 500, len(dates)-1])
def test_concatenate_with_overlap(overlap_factor):
    # test that lightcurves are being concatenated correctly both with and
    # without any overlap
    lc1 = sunpy.lightcurve.LightCurve.create(base_input, index=dates)
    dt = dates[1] - dates[0]
    # create a new lc that is shifted in time so there is some overlap
    lc2 = sunpy.lightcurve.LightCurve.create(base_input,
                                             index=[t + (dates[-1] - dates[0]) - (overlap_factor-1) * dt for t in dates])
    concat_lc = lc1.concatenate(lc2)
    assert len(concat_lc.data) == (len(lc1.data) + len(lc2.data) - overlap_factor)
    # check that the times are correct
    assert np.all(concat_lc.data.index[0:len(dates)] == lc1.data.index)
    # check that the original data is still there
    assert np.all(concat_lc.data.index[-len(lc2.data)+overlap_factor:] == lc2.data.index[overlap_factor:])
    # check that the new data is there
    assert np.all(concat_lc.data[-len(lc2.data)+overlap_factor:] == lc2.data[overlap_factor:])


def test_concatenate_meta():
    # check that meta data is also being added.
    eve = sunpy.lightcurve.EVELightCurve.create(EVE_AVERAGES_CSV)

    lc1 = sunpy.lightcurve.LightCurve.create(base_input, index=dates)
    dt = dates[1] - dates[0]
    # create a new lc that is shifted in time so there is some overlap
    lc2 = sunpy.lightcurve.LightCurve.create(base_input,
                                             index=[t + (dates[-1] - dates[0]) - 1 * dt for t in
                                                    dates])
    concat_lc = lc1.concatenate(lc2)
    assert len(concat_lc.meta) == len(lc1.meta) + len(lc2.meta)


def test_concatenate_fail():
    # check that concatenate throws an error when trying to concatenate
    # two different lightcurve classes
    lc1 = sunpy.lightcurve.GOESLightCurve.create(base_input, index=dates)
    eve = sunpy.lightcurve.EVELightCurve.create(EVE_AVERAGES_CSV)
    with pytest.raises(TypeError):
        lc1.concatenate(eve)


