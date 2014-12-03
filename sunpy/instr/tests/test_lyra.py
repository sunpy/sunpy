from __future__ import absolute_import
from __future__ import division

import tempfile
import os.path
import pytest
import datetime

import numpy as np

from sunpy.time import parse_time
from sunpy import config
from sunpy.instr import lyra

# Define location for test LYTAF database files
TEST_DATA_PATH = os.path.join(config.get("downloads", "download_dir"), "test")
# Define some test data for test_remove_lyra_artifacts()
TIME = np.array([datetime.datetime(2013, 2, 1)+datetime.timedelta(minutes=i)
                 for i in range(120)])
CHANNELS = [np.zeros(len(TIME))+0.4, np.array(len(TIME))+0.1]
EMPTY_LYTAF = np.empty((0,), dtype=[("insertion_time", object),
                                    ("begin_time", object),
                                    ("reference_time", object),
                                    ("end_time", object),
                                    ("event_type", object),
                                    ("event_definition", object)])
LYTAF = np.append(EMPTY_LYTAF,
                  np.array([(datetime.datetime.utcfromtimestamp(1371459961),
                             datetime.datetime.utcfromtimestamp(1359677220),
                             datetime.datetime.utcfromtimestamp(1359677250),
                             datetime.datetime.utcfromtimestamp(1359677400),
                             "LAR", "Large Angle Rotation.")],
                           dtype=EMPTY_LYTAF.dtype))
LYTAF = np.append(LYTAF,
                  np.array([(datetime.datetime.utcfromtimestamp(1371460063),
                             datetime.datetime.utcfromtimestamp(1359681764),
                             datetime.datetime.utcfromtimestamp(1359682450),
                             datetime.datetime.utcfromtimestamp(1359683136),
                             "UV occ.", "Occultation in the UV spectrum.")],
                           dtype=LYTAF.dtype))

@pytest.mark.online
def test_split_series_using_lytaf():
    '''test the downloading of the LYTAF file and subsequent queries'''
    tmp_dir = tempfile.mkdtemp()
    lyra.download_lytaf_database(lytaf_dir=tmp_dir)
    assert os.path.exists(os.path.join(tmp_dir, 'annotation_ppt.db'))

    #test split_series_using_lytaf
    #construct a dummy signal for testing purposes
    basetime = parse_time('2010-06-13 02:00')
    seconds = 3600
    dummy_time = [basetime + datetime.timedelta(0, s) for s in range(seconds)]
    dummy_data = np.random.random(seconds)

    lytaf = lyra.extract_lytaf_events('2010-06-13 02:00', '2010-06-13 06:00',
                                      lytaf_path=tmp_dir,
                                      combine_files=["ppt"])
    split = lyra.split_series_using_lytaf(dummy_time, dummy_data, lytaf)
    assert type(split) == list
    assert len(split) == 4
    assert split[0]['subtimes'][0] == datetime.datetime(2010, 6, 13, 2, 0)
    assert split[0]['subtimes'][-1] == datetime.datetime(2010, 6, 13, 2, 7, 2)
    assert split[3]['subtimes'][0] == datetime.datetime(2010, 6, 13, 2, 59, 41)
    assert split[3]['subtimes'][-1] == datetime.datetime(2010, 6, 13, 2, 59, 58)

def test_remove_lyra_artifacts_1():
    """Test remove_lyra_artifacts() with default values."""
    # Run remove_lyra_artifacts
    time_test, channels_test = lyra.remove_lyra_artifacts(
        TIME, channels=CHANNELS, lytaf_path=TEST_DATA_PATH,
        force_use_local_lytaf=True)
    # Generated expected result
    bad_indices = np.logical_or(
        np.logical_and(TIME >= LYTAF["begin_time"][0],
                       TIME <= LYTAF["end_time"][0]),
        np.logical_and(TIME >= LYTAF["begin_time"][1],
                       TIME <= LYTAF["end_time"][1]))
    bad_indices = np.arange(len(TIME))[bad_indices]
    time_expected = np.delete(TIME, bad_indices)
    channels_expected = [np.delete(CHANNELS[0], bad_indices),
                         np.delete(CHANNELS[1], bad_indices)]
    # Assert test values are same as expected
    assert time_test.all() == time_expected.all()
    assert (channels_test[0]).all() == (channels_expected[0]).all()
    assert (channels_test[1]).all() == (channels_expected[1]).all()

def test_remove_lyra_artifacts_2():
    """Test remove_lyra_artifacts() with user artifacts found and not found."""
    # Run remove_lyra_artifacts
    time_test, channels_test, artifacts_status_test = \
      lyra.remove_lyra_artifacts(
          TIME, channels=CHANNELS, artifacts=["LAR", "Offpoint"],
          return_artifacts=True, lytaf_path=TEST_DATA_PATH,
          force_use_local_lytaf=True)
    # Generated expected result
    bad_indices = np.logical_and(TIME >= LYTAF["begin_time"][0],
                                 TIME <= LYTAF["end_time"][0])
    bad_indices = np.arange(len(TIME))[bad_indices]
    time_expected = np.delete(TIME, bad_indices)
    channels_expected = [np.delete(CHANNELS[0], bad_indices),
                         np.delete(CHANNELS[1], bad_indices)]
    artifacts_status_expected = {"lytaf": LYTAF, "removed": LYTAF[0],
                                 "not_removed": LYTAF[1],
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

def test_remove_lyra_artifacts_3():
    """Test remove_lyra_artifacts() with no user artifacts found."""
    # Run remove_lyra_artifacts
    time_test, channels_test, artifacts_status_test = \
      lyra.remove_lyra_artifacts(
          TIME, channels=CHANNELS, artifacts=["Offpoint"],
          return_artifacts=True, lytaf_path=TEST_DATA_PATH,
          force_use_local_lytaf=True)
    # Generated expected result
    time_expected = TIME
    channels_expected = CHANNELS
    artifacts_status_expected = {"lytaf": LYTAF, "removed": EMPTY_LYTAF,
                                 "not_removed": LYTAF,
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

def test_remove_lyra_artifacts_4():
    """Test if correct errors are raised by remove_lyra_artifacts()."""
    with pytest.raises(TypeError):
        lyra.remove_lyra_artifacts(TIME, artifacts=[6],
                                   lytaf_path=TEST_DATA_PATH,
                                   force_use_local_lytaf=True)
    with pytest.raises(TypeError):
        lyra.remove_lyra_artifacts(TIME, channels=6,
                                   lytaf_path=TEST_DATA_PATH,
                                   force_use_local_lytaf=True)

def test_extract_lytaf_events():
    """Test if LYTAF events are correctly downloaded and read in."""
    # Run extract_combined_lytaf
    lytaf_test = lyra.extract_lytaf_events("2008-01-01", "2014-01-01",
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
        lytaf_test = lyra.extract_lytaf_events("2008-01-01", "2014-01-01",
                                               lytaf_path="test_data",
                                               combine_files=["gigo"],
                                               force_use_local_lytaf=True)
