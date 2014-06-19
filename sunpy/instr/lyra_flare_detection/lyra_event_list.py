# Outline of code
#-------Prepare data-------
# 1.  Create LYRALightcurve object for given time interval from level 3
#     fits files (1-min averaged).
# 2.  Download LYTAF file from server.  Add keyword to allow user to
#     define local file.
# 3.  Remove data during UV occulations, offpoints, calibrations, and
#     large angle rotations.  (Also remove recovery times?)
# 4.  Download Ingolf's comb-like masks.  Add keyword so they can be
#     defined locally.
# 5.  Use the masks to remove comb-like features.
#------Find possible start times flares--------
# 6.  Take derivative of time series.
# 7.  Take where derivative is positive.
# 8.  Take where derivative is negative.
# 9.  Copy time array and shift by 4 elements.
# 10. Subtract one from other and find where result = 4 minutes.
# 11. Define earlier time corresponding to previous step as preliminary
#     start times.
#-----Find end and peak times---------
# 12. Loop through start times.
# 13. Find next time where derivative in negative.
# 14. While end time criterion not reached, loop through following times
#     where derivative is negative and find when end criterion is
#     reached.
# 15. Define peak time as time of maximum between start and end time.
#-----Look for flares during decay phase------
# 16. If there are possible start times between peak and end time,
#     repeat steps 13-16.
from __future__ import absolute_import
from __future__ import division

import os.path
from datetime import datetime

import numpy as np
import sqlite3
from itertools import chain
from sunpy.time import parse_time

RISE_FACTOR = 1.4
FALL_FACTOR = 0.5

def extract_combined_lytaf(tstart, tend,
                           lytaf_path=os.path.join(os.path.curdir, "data"),
                           combine_files=["lyra", "manual", "ppt", "science"]):
    
    """
    Extracts combined lytaf file for given time range.

    Given a time range defined by start_time and end_time, this
    function extracts the segments of each LYRA annotation file and
    combines them.

    Parameters
    ----------
    start_time : datetime object or string
                 Start time of period for which annotation file is
                 required.
    end_time : datetime object or string
               End time of period for which annotation file is
               required.
    combine_files : (optional) list of strings
                    States which LYRA annotation files are to be
                    combined.
                    Default is all four, i.e. lyra, manual, ppt,
                    science.
                    See Notes section for an explanation of each.

    Returns
    -------
    lytaf : numpy record array containing the various parameters stored
            in the LYTAF files.

    Notes
    -----
    There are four LYRA annotation files which mark different types of
    events or artifacts in the data.  They are named
    annotation_suffix.db where suffix is a variable equalling either
    lyra, manual, ppt, or science.
    annotation_lyra.db : contains entries regarding possible effects to
                        the data due to normal operation of LYRA
                        instrument.
    annotation_manual.db : contains entries regarding possible effects
                           to the data due to unusual or manually
                           logged events.
    annotation_ppt.db : contains entries regarding possible effects to
                        the data due to pointing or positioning of
                        PROBA2.
    annotation_science.db : contains events in the data scientifically
                            interesting, e.g. flares.

    References
    ----------
    Further documentation: http://proba2.oma.be/data/TARDIS

    Examples
    --------
    
    """
    # Check inputs
    # Check start_time is a date string or datetime object
    if type(tstart) is str:
        tstart = parse_time(tstart)
    if type(tstart) is not datetime:
        raise TypeError("tstart must be a date string or datetime object")
    # Check start_time is a date string or datetime object
    if type(tend) is str:
        tend = parse_time(tend)
    if type(tend) is not datetime:
        raise TypeError("tend must be a date string or datetime object")
    # Check combine_files contains correct inputs
    if not all(suffix in ["lyra", "manual", "ppt", "science"]
               for suffix in combine_files):
        raise TypeError("Elements in combine_files must be strings equalling "
                        "'lyra', 'manual', 'ppt', or 'science'.")
    # Remove any duplicates from combine_files input
    combine_files = list(set(combine_files))
    combine_files.sort()
    # Convert input times to UNIX timestamp format since this is the
    # time format in the annotation files
    tstart_uts = tstart.strftime("%s")
    tend_uts = tend.strftime("%s")

    # Define numpy record array which will hold the information from
    # the annotation file.
    lytaf = np.empty((0,), dtype=[("insertion_time", object),
                               ("begin_time", object),
                               ("reference_time", object),
                               ("end_time", object),
                               ("event_type", object),
                               ("event_definition", object)])
    # Access annotation files
    for i, suffix in enumerate(combine_files):
        # Open SQLITE3 annotation files
        connection = sqlite3.connect(
            os.path.join(lytaf_path, "annotation_{0}.db".format(suffix)))
        # Create cursor to manipulate data in annotation file
        cursor = connection.cursor()
        # Select and extract the data from event table within file within
        # given time range
        cursor.execute("select insertion_time, begin_time, reference_time, "
                       "end_time, eventType_id from event where end_time >= "
                       "{0} and begin_time <= {1}".format(tstart_uts, tend_uts))
        event_rows = cursor.fetchall()
        # Select and extract the event types from eventType table
        cursor.row_factory = sqlite3.Row
        cursor.execute("select * from eventType")
        eventType_rows = cursor.fetchall()
        eventType_id = []
        eventType_type = []
        eventType_definition = []
        for eventType_row in eventType_rows:
            eventType_id.append(eventType_row["id"])
            eventType_type.append(eventType_row["type"])
            eventType_definition.append(eventType_row["definition"])
        # Enter desired information into the lytaf numpy record array
        for event_row in event_rows:
            id_index = eventType_id.index(event_row[4])
            lytaf = np.append(lytaf,
                              np.array((datetime.fromtimestamp(event_row[0]),
                                        datetime.fromtimestamp(event_row[1]),
                                        datetime.fromtimestamp(event_row[2]),
                                        datetime.fromtimestamp(event_row[3]),
                                        eventType_type[id_index],
                                        eventType_definition[id_index]),
                                        dtype=lytaf.dtype))
        # Close file
        cursor.close()
        connection.close()
    # Delete initial empty entry in lytaf
    np.delete(lytaf, 0)
    # Sort lytaf in ascending order of begin time
    np.recarray.sort(lytaf, order="begin_time")    

    #return event_rows, eventType_rows
    return lytaf
    
def find_lyra_events(flux, time):
    """
    Finds events in a times series satisfying LYRA event definitions.

    This function finds events in an input time series which satisfy
    the LYRA event definitions and returns the start, peak and end
    times.  For LYRA event definitions, see Notes section of this
    docstring.

    Parameters
    ----------
    flux : ndarray/array-like convertible to float64, e.g. np.array, list
           Contains flux/irradiance measurements
    time : ndarray/array-like of of datetime objects, e.g. np.array, list
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
    # Ensure inputs are of correct type
    flux = np.asanyarray(flux, dtype="float64")
    if (np.array([type(t) for t in time]) == datetime).all():
        time = np.asanyarray(time)
    else:
        # If elements of time are not datetime objects, try converting.
        try:
            time = np.array([datetime(t) for t in time])
        except TypeError:
            try:
                # If cannot be converted simply, elements may be strings
                # Try converting to datetime using sunpy.time.parse_time
                time = np.array([parse_time(t) for t in time])
            except:
                # Otherwise raise error telling user to input an array
                # of datetime objects.
                raise IOError("time must be an array or array-like of datetime"
                              " objects or valid time strings.")
        else:
            raise IOError("time must be an array or array-like of datetime"
                          " objects or valid time strings.")
    # Define variables to be used later
    flare_indices = []
    flare_times = []
    # Get LYTAF file for given time range
    #lytaf = extract_combined_lytaf(time[0], time[-1])
    # Find events in lytaf which are to be removed from time series.
    #artifacts = np.logical_or(lytaf["event_type"] == u'UV occ.',
    #                          lytaf["event_type"] == u'Offpoint',
    #                          lytaf["event_type"] == u'LAR',
    #                          lytaf["event_type"] == u'Calibration')
    # Remove periods corresponding to artifacts from flux and time arrays
    #for artifact in np.arange(len(artifacts))[artifacts]:
    #    bad_period = np.logical_and(flux > lytaf["begin_time"][index],
    #                                flux < lytaf["end_time"][index])
    #    flux = np.delete(flux, bad_period)
    #    time = np.delete(time, bad_period)
    # get derivative of flux wrt time
    time_timedelta = time[1:-1]-time[0:-2]
    dt = np.zeros(len(time_timedelta))
    for i, t, in enumerate(time_timedelta):
        dt[i] = t.total_seconds()
    dfdt = np.gradient(flux[0:-2], dt)
    # Get locations where derivative is positive
    pos_deriv = np.where(dfdt > 0)[0]
    neg_deriv = np.where(dfdt < 0)[0]
    # Find difference between each time point and the one 4
    # observations ahead.
    time_timedelta4 = time[4:-1]-time[0:-5]
    dt4 = np.zeros(len(time_timedelta4))
    for i, t, in enumerate(time_timedelta4):
        dt4[i] = t.total_seconds()
    # Find all possible flare start times.
    i=0
    while i < len(pos_deriv)-3:
        print i
        # Start time criteria
        if (pos_deriv[i:i+4]-pos_deriv[i] == np.arange(4)).all() and \
          dt4[i] > 210 and dt4[i] < 270 and \
          flux[pos_deriv[i+3]]/flux[pos_deriv[i]] >= RISE_FACTOR:
            # Next, find index of flare end time.
            jj = np.where(neg_deriv > pos_deriv[i])[0]
            j = neg_deriv[jj[0]]
            end_index = np.where(flux[j:] <= \
                                 max(flux[pos_deriv[i]:j])- \
                                 (max(flux[pos_deriv[i]:j])-flux[pos_deriv[i]]) \
                                 *FALL_FACTOR)[0][0]+j
            # find index of peak time
            peak_index = np.where(
                flux == max(flux[pos_deriv[i]:end_index]))[0][0]
            # Record flare start, peak and end times
            flare_indices.append((pos_deriv[i], peak_index, end_index))
            flare_times.append((time[pos_deriv[i]],
                               time[peak_index], time[end_index]))
            # If the most recently found flare is during the decay phase
            # of another reset end time of previous flare to start time
            # of this flare.
            #if flare_indices[-2][2] > flare_indices[-1][0]:
            #    new_peak_index = np.where(flux == max(
            #        flux[flare_indices[-2][0]:flare_indices[-1][0]]))[0][0]
            #    flare_indices[-2] = (flare_indices[-2][0],
            #                         new_peak_index, flare_indices[-1][0])
            #    flare_times[-2] = (time[flare_indices[-2][0]],
            #                       time[new_peak_index],
            #                       time[flare_indices[-1][0]])
            # Finally, set principle iterator, i, to the peak of the
            # flare just found so that algorithm will start looking for
            # flares during the decay phase of this flare and beyond,
            # thereby skipping the rise phase of this flare.  This
            # ensures that flares during the decay phase are also
            # located.
            i = peak_index
        else:
            i = i+1

    return flare_times
