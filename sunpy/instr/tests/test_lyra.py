from __future__ import absolute_import, division

import tempfile
import os.path
import pytest
import datetime

import numpy as np
import pandas

from sunpy.data.test import rootdir
from sunpy.time import parse_time
from sunpy import lightcurve
from sunpy import timeseries
from sunpy.instr import lyra

from sunpy.extern.six.moves import range

# Define location for test LYTAF database files
TEST_DATA_PATH = rootdir

# Define some test data for test_remove_lytaf_events()
TIME = np.array([datetime.datetime(2013, 2, 1) + datetime.timedelta(minutes=i)
                 for i in range(120)])
CHANNELS = [np.zeros(len(TIME)) + 0.4, np.zeros(len(TIME)) + 0.1]
EMPTY_LYTAF = np.empty((0,), dtype=[("insertion_time", object),
                                    ("begin_time", object),
                                    ("reference_time", object),
                                    ("end_time", object),
                                    ("event_type", object),
                                    ("event_definition", object)])
LYTAF_TEST = np.append(
    EMPTY_LYTAF,
    np.array([(datetime.datetime.utcfromtimestamp(1371459961),
               datetime.datetime.utcfromtimestamp(1359677220),
               datetime.datetime.utcfromtimestamp(1359677250),
               datetime.datetime.utcfromtimestamp(1359677400),
               "LAR", "Large Angle Rotation.")],
             dtype=EMPTY_LYTAF.dtype))
LYTAF_TEST = np.append(
    LYTAF_TEST,
    np.array([(datetime.datetime.utcfromtimestamp(1371460063),
               datetime.datetime.utcfromtimestamp(1359681764),
               datetime.datetime.utcfromtimestamp(1359682450),
               datetime.datetime.utcfromtimestamp(1359683136),
               "UV occ.", "Occultation in the UV spectrum.")],
             dtype=LYTAF_TEST.dtype))


@pytest.mark.remote_data
def test_split_series_using_lytaf():
    '''test the downloading of the LYTAF file and subsequent queries'''
    tmp_dir = tempfile.mkdtemp()
    lyra.download_lytaf_database(lytaf_dir=tmp_dir)
    assert os.path.exists(os.path.join(tmp_dir, 'annotation_ppt.db'))

    # test split_series_using_lytaf
    # construct a dummy signal for testing purposes
    basetime = parse_time('2010-06-13 02:00')
    seconds = 3600
    dummy_time = [basetime + datetime.timedelta(0, s) for s in range(seconds)]
    dummy_data = np.random.random(seconds)

    lytaf_tmp = lyra.get_lytaf_events('2010-06-13 02:00', '2010-06-13 06:00',
                                      lytaf_path=tmp_dir,
                                      combine_files=["ppt"])
    split = lyra.split_series_using_lytaf(dummy_time, dummy_data, lytaf_tmp)
    assert type(split) == list
    assert len(split) == 4
    assert split[0]['subtimes'][0] == datetime.datetime(2010, 6, 13, 2, 0)
    assert split[0]['subtimes'][-1] == datetime.datetime(2010, 6, 13, 2, 7, 2)
    assert split[3]['subtimes'][0] == datetime.datetime(2010, 6, 13, 2, 59, 41)
    assert split[3]['subtimes'][-1] == datetime.datetime(2010, 6, 13, 2, 59, 58)

    # Test case when no LYTAF events found in time series.
    split_no_lytaf = lyra.split_series_using_lytaf(dummy_time,
                                                   dummy_data, LYTAF_TEST)
    assert type(split_no_lytaf) == list
    assert type(split_no_lytaf[0]) == dict
    assert not set(split_no_lytaf[0].keys()).symmetric_difference({'subtimes', 'subdata'})
    assert split_no_lytaf[0]["subtimes"] == dummy_time
    assert split_no_lytaf[0]["subdata"].all() == dummy_data.all()


@pytest.fixture
def get_lyradata(dtype):
    if dtype == 'lc':
        # Create sample LYRALightCurve
        lyrats = lightcurve.LYRALightCurve.create("2014-01-01")
        lyrats.data = pandas.DataFrame(index=TIME,
                                       data={"CHANNEL1": CHANNELS[0],
                                             "CHANNEL2": CHANNELS[1],
                                             "CHANNEL3": CHANNELS[0],
                                             "CHANNEL4": CHANNELS[1]})
    else:
        # Create sample TimeSeries
        lyrats = timeseries.TimeSeries(
            os.path.join(rootdir, 'lyra_20150101-000000_lev3_std_truncated.fits.gz'),
            source='LYRA')
        lyrats.data = pandas.DataFrame(index=TIME,
                                       data={"CHANNEL1": CHANNELS[0],
                                             "CHANNEL2": CHANNELS[1],
                                             "CHANNEL3": CHANNELS[0],
                                             "CHANNEL4": CHANNELS[1]})
    return lyrats


@pytest.mark.remote_data
@pytest.mark.parametrize('dtype', ['ts', 'lc'])
def test_remove_lytaf_events_from_timeseries(dtype):
    """Test if artefacts are correctly removed from a TimeSeries.
    Also test LYRALightCurve for backwards compatibility."""
    # Check correct errors are raised due to bad input
    with pytest.raises(AttributeError):
        ts_test = lyra.remove_lytaf_events_from_timeseries(
            [], lytaf_path=TEST_DATA_PATH, force_use_local_lytaf=True)

    # Run remove_artifacts_from_timeseries, returning artifact
    # status
    ts = get_lyradata(dtype)
    ts_test, artifact_status_test = \
        lyra.remove_lytaf_events_from_timeseries(
            ts, artifacts=["LAR", "Offpoint"], return_artifacts=True,
            lytaf_path=TEST_DATA_PATH, force_use_local_lytaf=True)
    # Generate expected data by calling _remove_lytaf_events and
    # constructing expected dataframe manually.
    time, channels, artifact_status_expected = lyra._remove_lytaf_events(
        ts.data.index, channels=[np.asanyarray(ts.data["CHANNEL1"]),
                                 np.asanyarray(ts.data["CHANNEL2"]),
                                 np.asanyarray(ts.data["CHANNEL3"]),
                                 np.asanyarray(ts.data["CHANNEL4"])],
        artifacts=["LAR", "Offpoint"], return_artifacts=True,
        lytaf_path=TEST_DATA_PATH, force_use_local_lytaf=True)
    dataframe_expected = pandas.DataFrame(index=time,
                                          data={"CHANNEL1": channels[0],
                                                "CHANNEL2": channels[1],
                                                "CHANNEL3": channels[2],
                                                "CHANNEL4": channels[3]})
    # Assert expected result is returned
    pandas.util.testing.assert_frame_equal(ts_test.data, dataframe_expected)
    assert artifact_status_test.keys() == artifact_status_expected.keys()
    np.testing.assert_array_equal(artifact_status_test["lytaf"],
                                  artifact_status_expected["lytaf"])
    np.testing.assert_array_equal(artifact_status_test["removed"],
                                  artifact_status_expected["removed"])
    np.testing.assert_array_equal(artifact_status_test["not_removed"],
                                  artifact_status_expected["not_removed"])
    assert artifact_status_test["not_found"] == \
        artifact_status_expected["not_found"]

    # Run remove_artifacts_from_timeseries, without returning
    # artifact status
    ts_test = \
        lyra.remove_lytaf_events_from_timeseries(
            ts, artifacts=["LAR", "Offpoint"],
            lytaf_path=TEST_DATA_PATH, force_use_local_lytaf=True)
    # Assert expected result is returned
    pandas.util.testing.assert_frame_equal(ts_test.data, dataframe_expected)


def test_remove_lytaf_events_1():
    """Test _remove_lytaf_events() with some artifacts found and others not."""
    # Run _remove_lytaf_events
    time_test, channels_test, artifacts_status_test = \
        lyra._remove_lytaf_events(
            TIME, channels=CHANNELS, artifacts=["LAR", "Offpoint"],
            return_artifacts=True, lytaf_path=TEST_DATA_PATH,
            force_use_local_lytaf=True)
    # Generated expected result
    bad_indices = np.logical_and(TIME >= LYTAF_TEST["begin_time"][0],
                                 TIME <= LYTAF_TEST["end_time"][0])
    bad_indices = np.arange(len(TIME))[bad_indices]
    time_expected = np.delete(TIME, bad_indices)
    channels_expected = [np.delete(CHANNELS[0], bad_indices),
                         np.delete(CHANNELS[1], bad_indices)]
    artifacts_status_expected = {"lytaf": LYTAF_TEST, "removed": LYTAF_TEST[0],
                                 "not_removed": LYTAF_TEST[1],
                                 "not_found": ["Offpoint"]}
    # Assert test values are same as expected
    assert time_test.all() == time_expected.all()
    assert (channels_test[0]).all() == (channels_expected[0]).all()
    assert (channels_test[1]).all() == (channels_expected[1]).all()
    assert artifacts_status_test.keys() == artifacts_status_expected.keys()
    np.testing.assert_array_equal(artifacts_status_test["lytaf"],
                                  artifacts_status_expected["lytaf"])
    np.testing.assert_array_equal(artifacts_status_test["removed"],
                                  artifacts_status_expected["removed"])
    np.testing.assert_array_equal(artifacts_status_test["not_removed"],
                                  artifacts_status_expected["not_removed"])
    assert artifacts_status_test["not_found"] == \
        artifacts_status_expected["not_found"]

    # Test that correct values are returned when channels kwarg not
    # supplied.
    # Run _remove_lytaf_events
    time_test, artifacts_status_test = \
        lyra._remove_lytaf_events(
            TIME, artifacts=["LAR", "Offpoint"],
            return_artifacts=True, lytaf_path=TEST_DATA_PATH,
            force_use_local_lytaf=True)
    # Assert test values are same as expected
    assert time_test.all() == time_expected.all()
    assert artifacts_status_test.keys() == artifacts_status_expected.keys()
    np.testing.assert_array_equal(artifacts_status_test["lytaf"],
                                  artifacts_status_expected["lytaf"])
    np.testing.assert_array_equal(artifacts_status_test["removed"],
                                  artifacts_status_expected["removed"])
    np.testing.assert_array_equal(artifacts_status_test["not_removed"],
                                  artifacts_status_expected["not_removed"])
    assert artifacts_status_test["not_found"] == \
        artifacts_status_expected["not_found"]


def test_remove_lytaf_events_2():
    """Test _remove_lytaf_events() with no user artifacts found."""
    # Run _remove_lytaf_events
    time_test, channels_test, artifacts_status_test = \
        lyra._remove_lytaf_events(
            TIME, channels=CHANNELS, artifacts="Offpoint",
            return_artifacts=True, lytaf_path=TEST_DATA_PATH,
            force_use_local_lytaf=True)
    # Generated expected result
    time_expected = TIME
    channels_expected = CHANNELS
    artifacts_status_expected = {"lytaf": LYTAF_TEST, "removed": EMPTY_LYTAF,
                                 "not_removed": LYTAF_TEST,
                                 "not_found": ["Offpoint"]}
    # Assert test values are same as expected
    assert time_test.all() == time_expected.all()
    assert (channels_test[0]).all() == (channels_expected[0]).all()
    assert (channels_test[1]).all() == (channels_expected[1]).all()
    assert artifacts_status_test.keys() == artifacts_status_expected.keys()
    np.testing.assert_array_equal(artifacts_status_test["lytaf"],
                                  artifacts_status_expected["lytaf"])
    np.testing.assert_array_equal(artifacts_status_test["removed"],
                                  artifacts_status_expected["removed"])
    np.testing.assert_array_equal(artifacts_status_test["not_removed"],
                                  artifacts_status_expected["not_removed"])
    assert artifacts_status_test["not_found"] == \
        artifacts_status_expected["not_found"]

    # Test correct values are returned when return_artifacts kwarg not
    # supplied.
    # Case 1: channels kwarg is True
    # Run _remove_lytaf_events
    time_test, channels_test = lyra._remove_lytaf_events(
        TIME, channels=CHANNELS, artifacts=["Offpoint"],
        lytaf_path=TEST_DATA_PATH, force_use_local_lytaf=True)
    assert time_test.all() == time_expected.all()
    assert (channels_test[0]).all() == (channels_expected[0]).all()
    assert (channels_test[1]).all() == (channels_expected[1]).all()
    # Case 2: channels kwarg is False
    # Run _remove_lytaf_events
    time_test = lyra._remove_lytaf_events(
        TIME, artifacts=["Offpoint"],
        lytaf_path=TEST_DATA_PATH, force_use_local_lytaf=True)
    assert time_test.all() == time_expected.all()


def test_remove_lytaf_events_3():
    """Test if correct errors are raised by _remove_lytaf_events()."""
    with pytest.raises(TypeError):
        lyra._remove_lytaf_events(TIME, channels=6, artifacts=["LAR"],
                                  lytaf_path=TEST_DATA_PATH,
                                  force_use_local_lytaf=True)
    with pytest.raises(ValueError):
        lyra._remove_lytaf_events(TIME,
                                  lytaf_path=TEST_DATA_PATH,
                                  force_use_local_lytaf=True)
    with pytest.raises(TypeError):
        lyra._remove_lytaf_events(TIME, artifacts=[6],
                                  lytaf_path=TEST_DATA_PATH,
                                  force_use_local_lytaf=True)
    with pytest.raises(ValueError):
        lyra._remove_lytaf_events(TIME,
                                  artifacts=["LAR", "incorrect artifact type"],
                                  lytaf_path=TEST_DATA_PATH,
                                  force_use_local_lytaf=True)


def test_get_lytaf_events():
    """Test if LYTAF events are correctly downloaded and read in."""
    # Run get_lytaf_events
    lytaf_test = lyra.get_lytaf_events("2008-01-01", "2014-01-01",
                                       lytaf_path=TEST_DATA_PATH,
                                       force_use_local_lytaf=True)
    # Form expected result of extract_combined_lytaf
    insertion_time = [datetime.datetime.utcfromtimestamp(1371459961),
                      datetime.datetime.utcfromtimestamp(1371460063),
                      datetime.datetime.utcfromtimestamp(1371460411),
                      datetime.datetime.utcfromtimestamp(1371460493),
                      datetime.datetime.utcfromtimestamp(1371460403),
                      datetime.datetime.utcfromtimestamp(1371470988),
                      datetime.datetime.utcfromtimestamp(1371211791),
                      datetime.datetime.utcfromtimestamp(1371212303)]
    begin_time = [datetime.datetime.utcfromtimestamp(1359677220),
                  datetime.datetime.utcfromtimestamp(1359681764),
                  datetime.datetime.utcfromtimestamp(1360748513),
                  datetime.datetime.utcfromtimestamp(1361115900),
                  datetime.datetime.utcfromtimestamp(1361980964),
                  datetime.datetime.utcfromtimestamp(1368581100),
                  datetime.datetime.utcfromtimestamp(1371032084),
                  datetime.datetime.utcfromtimestamp(1371158167)]
    reference_time = [datetime.datetime.utcfromtimestamp(1359677250),
                      datetime.datetime.utcfromtimestamp(1359682450),
                      datetime.datetime.utcfromtimestamp(1360751528),
                      datetime.datetime.utcfromtimestamp(1361116200),
                      datetime.datetime.utcfromtimestamp(1361983979),
                      datetime.datetime.utcfromtimestamp(1368582480),
                      datetime.datetime.utcfromtimestamp(1371045475),
                      datetime.datetime.utcfromtimestamp(1371162600)]
    end_time = [datetime.datetime.utcfromtimestamp(1359677400),
                datetime.datetime.utcfromtimestamp(1359683136),
                datetime.datetime.utcfromtimestamp(1360754543),
                datetime.datetime.utcfromtimestamp(1361116320),
                datetime.datetime.utcfromtimestamp(1361986994),
                datetime.datetime.utcfromtimestamp(1368583080),
                datetime.datetime.utcfromtimestamp(1371050025),
                datetime.datetime.utcfromtimestamp(1371167100)]
    event_type = [u"LAR", u"UV occ.", u"Vis LED on", u"M Flare", u"UV LED on",
                  u"X Flare", u"Off-limb event", u"Unexplained feature"]
    event_description = [u"Large Angle Rotation.",
                         u"Occultation in the UV spectrum.",
                         u"Visual LED is turned on.",
                         u"M class solar flare.",
                         u"UV LED is turned on.",
                         u"X class solar flare.",
                         u"Off-limb event in SWAP.",
                         u"Unexplained feature."]
    lytaf_expected = np.empty((8,), dtype=[("insertion_time", object),
                                           ("begin_time", object),
                                           ("reference_time", object),
                                           ("end_time", object),
                                           ("event_type", object),
                                           ("event_definition", object)])
    lytaf_expected["insertion_time"] = insertion_time
    lytaf_expected["begin_time"] = begin_time
    lytaf_expected["reference_time"] = reference_time
    lytaf_expected["end_time"] = end_time
    lytaf_expected["event_type"] = event_type
    lytaf_expected["event_definition"] = event_description
    # Assert that extract_combined_lytaf gives the right result
    np.testing.assert_array_equal(lytaf_test, lytaf_expected)

    # Check correct error is raised if names of different lytaf files
    # are incorrectly input.
    with pytest.raises(ValueError):
        lytaf_test = lyra.get_lytaf_events("2008-01-01", "2014-01-01",
                                           lytaf_path="test_data",
                                           combine_files=["gigo"],
                                           force_use_local_lytaf=True)


def test_get_lytaf_event_types():
    """Test that LYTAF event types are printed."""
    lyra.get_lytaf_event_types(lytaf_path=TEST_DATA_PATH)


def test_lytaf_event2string():
    """Test _lytaf_event2string() associates correct numbers and events."""
    out_test = lyra._lytaf_event2string(list(range(12)))
    assert out_test == ['LAR', 'N/A', 'UV occult.', 'Vis. occult.', 'Offpoint',
                        'SAA', 'Auroral zone', 'Moon in LYRA', 'Moon in SWAP',
                        'Venus in LYRA', 'Venus in SWAP']
    out_test_single = lyra._lytaf_event2string(1)
    assert out_test_single == ['LAR']


def test_prep_columns():
    """Test whether _prep_columns correctly prepares data."""
    # Generate simple input data
    time_input = TIME[0:2]
    channels_input = [CHANNELS[0][0:2], CHANNELS[1][0:2]]
    filecolumns_input = ["time", "channel0", "channel1"]

    # Test case when channels and filecolumns are supplied by user.
    string_time_test, filecolumns_test = lyra._prep_columns(
        time_input, channels_input, filecolumns_input)
    # Generate expected output and verify _prep_columns() works
    string_time_expected = np.array([t.strftime("%Y-%m-%dT%H:%M:%S.%f")
                                     for t in time_input])
    filecolumns_expected = ["time", "channel0", "channel1"]
    np.testing.assert_array_equal(string_time_test, string_time_expected)
    assert filecolumns_test == filecolumns_expected

    # Test case when channels supplied by user by not filecolumns
    string_time_test, filecolumns_test = lyra._prep_columns(time_input,
                                                            channels_input)
    np.testing.assert_array_equal(string_time_test, string_time_expected)
    assert filecolumns_test == filecolumns_expected

    # Test case when neither channels nor filecolumns supplied by user
    string_time_test, filecolumns_test = lyra._prep_columns(time_input)
    np.testing.assert_array_equal(string_time_test, string_time_expected)
    assert filecolumns_test == ["time"]

    # Test correct exceptions are raised
    with pytest.raises(TypeError):
        string_time_test, filecolumns_test = lyra._prep_columns(
            time_input, channels_input, ["channel0", 1])
    with pytest.raises(ValueError):
        string_time_test = lyra._prep_columns(time_input,
                                              filecolumns=filecolumns_input)
