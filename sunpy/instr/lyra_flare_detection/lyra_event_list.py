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
# 10.  Subtract one from other and find where result = 4 minutes.
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
from sunpy.time import parse_time

import sqlite3
from itertools import chain

def extract_combined_lytaf(tstart, tend,
                           lytaf_path=os.path.join(os.path.curdir, "data")
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
    lytaf : dictionary containing the various parameters stored in the
            LYTAF files.

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
    if type(tstart) is not datetime.datetime:
        raise TypeError("start_time must be a date string or datetime object")
    # Check start_time is a date string or datetime object
    if type(tend) is str:
        tend = parse_time(tend)
    if type(tend) is not datetime.datetime:
        raise TypeError("end_time must be a date string or datetime object")
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
    tstart_uts = tstart - datetime(1973,1,1).total_seconds()
    tend_uts = tend - datetime(1973,1,1).total_seconds()

    # Access annotation files
    # Define list to hold data from event tables in annotation files.
    # First element is defined as empty list for convenience in for loop below.
    # This will be deleted afterwards.
    event_rows = [[]]
    for i, suffix in enumerate(combine_files):
        # Open SQLITE3 annotation files
        connection = sqlite3.connect(os.path.join(lytaf_path, "annotation_"
                                                     + suffix + ".db"))
        # Create cursor to manipulate data in annotation file
        cursor = connection.cursor()
        # Select and extract the data from event table within file within
        # given time range
        cursor.execute("select insertion_time, begin_time, reference_time, "
                       "end_time, eventType_id from event where begin_time >= "
                       + tstart_uts + " and begin_time <= " + tend_uts)
        event_rows.append(cursor.fetchall())
        # Select and extract the event types from eventType table
        cursor.row_factory = sqlite3.Row
        cursor.execute("select * from eventType")
        eventType_rows = cursor.fetchall()
        eventType_id = []
        eventType_type = []
        eventType_definition = []
        for j, row in enumerate(eventType_rows):
            eventType_id.append(row["id"])
            eventType_type.append(row["type"])
            eventType_definition.append(row["definition"])
        # Replace event type IDs with event type description
        for j, row in enumerate(event_rows[i+1]):
            _id = eventType_id.index(event_rows[i+1][j][4])
            event_rows[i+1][j] = (event_rows[i+1][j][0], event_rows[i+1][j][1],
                                  event_rows[i+1][j][2], event_rows[i+1][j][3],
                                  eventType_rows[_id][1],
                                  eventType_rows[_id][2])
        # Close file
        cursor.close()
        connection.close()
    # Delete empty string as first element in event_rows
    del(event_rows[0])
    
    # Format and Output data
    # Make event_rows a unified list rather than a list of lists
    # where each sublist corresponds to an annotation file.
    event_rows = list(chain.from_iterable(event_rows))
    # Sort arrays in order of increasing start time.
    event_rows.sort(key=lambda tup: tup[1])

    # Convert times to datetime objects
    # INSERT CODE HERE

    # Create dictionary to hold results
    lytaf = {"insertion_time": [],
             "begin_time": [],
             "reference_time": [],
             "end_time": [],
             "event_type": [],
             "event_definition": []
            }
    # Put results into dictionary using for loop
    for i, row in enumerate(event_rows):
        lytaf["insertion_time"].append(datetime.fromtimestamp(
            event_rows[i][0]))
        lytaf["begin_time"].append(datetime.fromtimestamp(event_rows[i][1]))
        lytaf["reference_time"].append(datetime.fromtimestamp(
            event_rows[i][2]))
        lytaf["end_time"].append(datetime.fromtimestamp(event_rows[i][3]))
        lytaf["event_type"].append(event_rows[i][4])
        lytaf["event_definition"].append(event_rows[i][5])

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
    # Get LYTAF file for given time range
    lytaf = extract_combined_lytaf(time[0], time[-1])
    # Extract




#Get list of lists
#Replace eventType_id with eventtype for each tuple in each sublist using nested for loops.
#Unpack list of lists
#Order by event start time
#Generate output structure ?Dictionary?
#return
