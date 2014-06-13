from __future__ import absolute_import
from __future__ import division

from datetime import datetime
import pytest
import unittest

import lyra_event_list as lel

def test_extract_combined_lytaf():
    # Run extract_combined_lytaf
    lytaf_test = lel.extract_combined_lytaf("2008-01-01", "2014-01-01",
                                            lytaf_path="test/test_data")
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
    lytaf_expected = {"insertion_time": insertion_time,
                      "begin_time": begin_time,
                      "reference_time": reference_time,
                      "end_time": end_time,
                      "event_type": event_type,
                      "event_description": event_description}
    # Assert that extract_combined_lytaf gives the right result
    assert lytaf_test == lytaf_expected
    
