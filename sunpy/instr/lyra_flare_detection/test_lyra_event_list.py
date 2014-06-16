from __future__ import absolute_import
from __future__ import division

from datetime import datetime
import pytest
import unittest

import numpy as np
from numpy.testing import assert_array_equal

import lyra_event_list

def test_extract_combined_lytaf():
    # Run extract_combined_lytaf
    lytaf_test = lyra_event_list.extract_combined_lytaf("2008-01-01", "2014-01-01",
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
    event_type = [u'LAR', u'UV occ.', u'Vis LED on', u'M Flare', u'UV LED on',
                  u'X Flare', u'Off-limb event', u'Unexplained feature']
    event_definition = [u'Large Angle Rotation.',
                        u'Occultation in the UV spectrum.',
                        u'Visible LED is turned on.', 
                        u'M class solar flare.',
                        u'UV LED is turned on.',
                        u'X class solar flare.',
                        u'Off-limb event in SWAP.',
                        u'Unexplained feature.']
    n = len(begin_time)
    lytaf_expected = np.empty((n,), dtype=[("insertion_time", object),
                                        ("begin_time", object),
                                        ("reference_time", object),
                                        ("end_time", object),
                                        ("event_type", object),
                                        ("event_definition", object)])
    for i in range(n):
        lytaf_expected[i] = (insertion_time[i], begin_time[i],
                             reference_time[i], end_time[i],
                             event_type[i], event_definition[i])
    return lytaf_test, lytaf_expected
    # check function gives expected result
    assert lytaf_test.dtype.names == lytaf_expected.dtype.names
    assert (lytaf_test["insertion_time"] ==
            lytaf_expected["insertion_time"]).all()
    assert (lytaf_test["begin_time"] == lytaf_expected["begin_time"]).all()
    assert (lytaf_test["reference_time"] ==
            lytaf_expected["reference_time"]).all()
    assert (lytaf_test["end_time"] == lytaf_expected["end_time"]).all()
    assert (lytaf_test["event_type"] == lytaf_expected["event_type"]).all()
    assert (lytaf_test["event_definition"] ==
            lytaf_expected["event_definition"]).all()
