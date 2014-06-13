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

def find_lyra_events(flux, time):
    """
    Finds events in a times series satisfying LYRA event definitions.

    This function finds events in an input time series which satisfy
    the LYRA event definitions and returns the start, peak and end
    times.  For LYRA event definitions, see Notes section of this
    docstring.

    Parameters
    ----------
    flux : numpy float64 array
           Contains flux/irradiance measurements
    time : sunpy datetime object or string array/list
           Contains measurement times corresponding to each element in
           flux.  Must be same length as flux.

    Returns
    -------
    event_list :

    Notes
    -----
    Start time:
    1) There must be a continuous increase in 1-minute-averaged data
    over 4 minutes.
    2) The flux in the 4th minute must be at least 1.4 times the flux
    in the first minute.
    End time:
    1) The end time is when the flux falls to half-way between the peak
    and initial fluxes.

    References
    ---------
    http://www.swpc.noaa.gov/ftpdir/indices/events/README

    Examples
    --------                
    
    """
    # Define variables to be used later
    flare_indices = []
    flare_times = []
    # Get LYTAF file for given time range
    lytaf = extract_combined_lytaf(time[0], time[-1])
    # Find events in lytaf which are to be removed from time series.
    artifacts = np.logical_or(lytaf["event_type"] == u'UV occ.',
                              lytaf["event_type"] == u'Offpoint',
                              lytaf["event_type"] == u'LAR',
                              lytaf["event_type"] == u'Calibration')
    # Remove periods corresponding to artifacts from flux and time arrays
    for artifact in np.arange(len(artifacts))[artifacts]:
        bad_period = np.logical_and(flux > lytaf["begin_time"][index],
                                    flux < lytaf["end_time"][index])
        flux = np.delete(flux, bad_period)
        time = np.delete(time, bad_period)
    # get derivative of flux wrt time
    dt = time[0:-2] - time[1:-1]
    dfdt = np.gradient(flux[0:-2], dt)
    dfdt = np.append(dfdt, 0)
    # Get locations where derivative is positive
    pos_deriv = np.where(dfdt > 0)[0]
    neg_deriv = np.where(dfdt < 0)[0]
    # Find difference between each time point and the one 4
    # observations ahead.
    dt4 = time[0:-5] - time[4:-1]
    # Find all possible flare start times.
    for i in np.arange(len(pos_deriv)):
        # Start time criteria
        if pos_deriv[i:i+4]-pos_deriv[i] == np.arange(5) and dt4[i] == 4 and \
          flux[pos_deriv[i+4]] / flux[pos_deriv[i]] >= RISE_FACTOR:
            # Next, find index of flare end time.
            jj = np.where(neg_deriv > pos_deriv[i])[0]
            j = neg_deriv[jj[0]]
            end_index = np.where(flux[j:] <=
                                 max(flux[i:j]) - (max(flux[i:j])-flux[i]) *
                                 FALL_FACTOR)[0][0] + j
            # find index of peak time
            peak_index = np.where(flux == max(flux[pos_deriv[i]:j]))[0][0]
            # Record flare start, peak and end times
            flare_indices.append(pos_deriv[i]], peak_index, end_index)
            flare_times.append(time[pos_deriv[i]],
                               time[peak_index], time[end_index])
            # If the most recently found flare is during the decay phase
            # of another reset end time of previous flare to start time
            # of this flare.
            if flare_indices[-2][2] > flare_indices[-1][0]:
                new_peak_index = np.where(flux == max(
                    flux[flare_indices[-2][0]:flare_indices[-1][0]]))[0][0]
                flare_indices[-2] = (flare_indices[-2][0],
                                     new_peak_index, flare_indices[-1][0])
                flare_times[-2] = (time[flare_indices[-2][0]],
                                   time[new_peak_index],
                                   time[flare_indices[-1][0]])
            # Finally, set principle iterator, i, to the peak of the
            # flare just found so that algorithm will start looking for
            # flares during the decay phase of this flare and beyond,
            # thereby skipping the rise phase of this flare.  This
            # ensures that flares during the decay phase are also
            # located.
            i = peak_index

    return flare_times

    
