from __future__ import absolute_import

from datetime import datetime

import sunpy
from sunpy.roi import roi


def test_roi_instance():
    region = roi(times=['2012-06-20 05:00', '2012-06-20 07:00'],
                 description='dummy_roi')
    assert isinstance(region, sunpy.roi.roi)


def test_roi_empty_instance():
    region = roi()
    assert (region.start_time is None)
    assert (region.end_time is None)


def test_roi_times_str():
    region = roi(times='2012-06-20 05:00')
    expected_time = datetime(2012, 6, 20, 5, 0)
    assert (region.start_time == expected_time)
    assert (region.end_time == expected_time)


def test_roi_times_list_one_element():
    region = roi(times=['2012-06-20 05:00'])
    expected_time = datetime(2012, 6, 20, 5, 0)
    assert (region.start_time == expected_time)
    assert (region.end_time == expected_time)


def test_roi_times_list_two_elements():
    region = roi(times=['2012-06-20 05:00', '2012-06-20 07:00'])
    expected_start_time = datetime(2012, 6, 20, 5, 0)
    expected_end_time = datetime(2012, 6, 20, 7, 0)
    assert (region.start_time == expected_start_time)
    assert (region.end_time == expected_end_time)


def test_roi_times_list_more_that_2_elements():
    region = roi(times=['2012-06-20 05:00', '2012-06-20 07:00',
                        '2012-06-20 09:00'])
    assert (region.start_time is None)
    assert (region.end_time is None)


def test_roi_description():
    region = roi(description='foo')
    assert isinstance(region, sunpy.roi.roi)
    assert (region.description == 'foo')


def test_roi_source():
    region = roi(source='foo')
    assert isinstance(region, sunpy.roi.roi)
    assert (region.source == 'foo')


def test_roi_time_range():
    region = roi(times=['2012-06-20 05:00', '2012-06-20 07:00'],
                 description='dummy_roi')
    assert isinstance(region.time_range(), sunpy.time.TimeRange)
