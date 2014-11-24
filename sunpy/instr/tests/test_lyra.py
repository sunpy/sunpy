from __future__ import absolute_import

import tempfile
import os.path
import pytest
import datetime

import numpy as np

from sunpy.time import TimeRange, parse_time
from sunpy.instr import lyra
from sunpy import config

TEST_DATA_PATH = os.path.join(config.get("downloads", "download_dir"), "test")

@pytest.mark.online
def test_lytaf_utils():
    '''test the downloading of the LYTAF file and subsequent queries'''
    tmp_dir = tempfile.mkdtemp()
    lyra.download_lytaf_database(lytaf_dir=tmp_dir)
    assert os.path.exists(os.path.join(tmp_dir,'annotation_ppt.db'))

    #try doing a query on the temporary database
    lar = lyra.get_lytaf_events(
        TimeRange('2010-06-13 02:00','2010-06-13 06:00'), lytaf_dir=tmp_dir)
    assert type(lar) == list
    assert type(lar[0]) == dict
    assert type(lar[0]['start_time']) == datetime.datetime
    assert type(lar[0]['end_time']) == datetime.datetime
    assert type(lar[0]['roi_description']) == str
    assert type(lar[0]['event_type_description']) == str
    assert lar[0]['start_time'] == parse_time('2010-06-13 02:07:04')
    assert lar[0]['end_time'] == parse_time('2010-06-13 02:10:04')
    assert lar[0]['event_type_description'] == 'LAR'

    #test split_series_using_lytaf
    #construct a dummy signal for testing purposes
    basetime = parse_time('2010-06-13 02:00')
    seconds = 3600
    dummy_time = [basetime + datetime.timedelta(0, s) for s in range(seconds)]
    dummy_data = np.random.random(seconds)

    split = lyra.split_series_using_lytaf(dummy_time, dummy_data, lar)
    assert type(split) == list
    assert len(split) == 4
    assert split[0]['subtimes'][0] == datetime.datetime(2010, 6, 13, 2, 0)
    assert split[0]['subtimes'][-1] == datetime.datetime(2010, 6, 13, 2, 7, 2)
    assert split[3]['subtimes'][0] == datetime.datetime(2010, 6, 13, 2, 59, 41)
    assert split[3]['subtimes'][-1] == datetime.datetime(2010, 6, 13, 2, 59, 58)

def test_remove_lyra_artifacts():
    """Test if times during LYTAF events are removed from LYRA timeseries."""

def test_get_lytaf_events():
    """Test if LYTAF events are correctly downloaded and read in."""
    # Run extract_combined_lytaf
    lytaf_test = lyra.get_lytaf_events("2008-01-01", "2014-01-01",
                                       lytaf_path=TEST_DATA_PATH)
    # Form expected result of extract_combined_lytaf
    insertion_time = [datetime.fromtimestamp(1371459961),
                      datetime.fromtimestamp(1371460063),
                      datetime.fromtimestamp(1371460411),
                      datetime.fromtimestamp(1371460493),
                      datetime.fromtimestamp(1371460403),
                      datetime.fromtimestamp(1371470988),
                      datetime.fromtimestamp(1371211791),
                      datetime.fromtimestamp(1371212303)]
    begin_time = [datetime.fromtimestamp(1359677220),
                  datetime.fromtimestamp(1359681764),
                  datetime.fromtimestamp(1360748513),
                  datetime.fromtimestamp(1361115900),
                  datetime.fromtimestamp(1361980964),
                  datetime.fromtimestamp(1368581100),
                  datetime.fromtimestamp(1371032084),
                  datetime.fromtimestamp(1371158167)]
    reference_time = [datetime.fromtimestamp(1359677250),
                      datetime.fromtimestamp(1359682450),
                      datetime.fromtimestamp(1360751528),
                      datetime.fromtimestamp(1361116200),
                      datetime.fromtimestamp(1361983979),
                      datetime.fromtimestamp(1368582480),
                      datetime.fromtimestamp(1371045475),
                      datetime.fromtimestamp(1371162600)]
    end_time = [datetime.fromtimestamp(1359677400),
                datetime.fromtimestamp(1359683136),
                datetime.fromtimestamp(1360754543),
                datetime.fromtimestamp(1361116320),
                datetime.fromtimestamp(1361986994),
                datetime.fromtimestamp(1368583080),
                datetime.fromtimestamp(1371050025),
                datetime.fromtimestamp(1371167100)]
    event_type = ["LAR", "UV occ.", "Vis LED on", "M Flare", "UV LED on",
                  "X Flare", "Off-limb event", "Unexplained feature"]
    event_description = ["Large Angle Rotation.",
                         "Occultation in the UV spectrum.",
                         "Visible LED is turned on.", 
                         "M class solar flare.",
                         "UV LED is turned on.",
                         "X class solar flare.",
                         "Off-limb event in SWAP.",
                         "Unexplained feature."]
    lytaf_expected = lytaf = np.empty((8,), dtype=[("insertion_time", object),
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
    lytaf_expected["event_decription"] = event_description
    # Assert that extract_combined_lytaf gives the right result
    assert lytaf_test == lytaf_expected

    # Check correct error is raised if names of different lytaf files
    # are incorrectly input.
    with pytest.raises(ValueError):
        lytaf_test = lel.extract_combined_lytaf("2008-01-01", "2014-01-01",
                                                lytaf_path="test_data",
                                                combine_files=["gigo"])
