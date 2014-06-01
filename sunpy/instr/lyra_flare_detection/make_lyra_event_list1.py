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
import datetime

import sqlite3
import numpy as np
from sunpy import parse_time

# Define location of LYRA Annotation Files
LYTAF_PATH = os.path.join(os.path.curdir, "data")

def extract_combined_lytaf(tstart, tend,
                           combine_files=["lyra", "manual", "ppt", "science"])
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
    if type(start_time) is str:
        start_time = parse_time(start_time)
    if type(start_time) is not datetime.datetime:
        raise TypeError("start_time must be a date string or datetime object")
    # Check start_time is a date string or datetime object
    if type(end_time) is str:
        end_time = parse_time(end_time)
    if type(end_time) is not datetime.datetime:
        raise TypeError("end_time must be a date string or datetime object")
    # Check combine_files contains correct inputs
    if not all(suffix in ["lyra", "manual", "ppt", "science"]
               for suffix in combine_files):
        raise TypeError("Elements in combine_files must be strings equalling "
                        "'lyra', 'manual', 'ppt', or 'science'.")
    # Remove any duplicates from combine_files input
    combine_files = list(set(combine_files))

    # Access annotation files
    insertion_time = []
    begin_time = []
    reference_time = []
    end_time = []
    eventType_id = []
    _id = []
    _type = []
    _definition = []
    for suffix in combine_files:
        # Open SQLITE3 annotation files
        connection = sqlite3.connect(os.path.join(LYTAF_PATH, "annotation_"
                                                     + suffix + ".db"))
        # Create cursor to manipulate data in annotation file
        cursor = connection.cursor()
        # Select and extract the data from event table within file
        cursor.row_factory = sqlite3.Row
        cursor.execute("SELECT * FROM event")
        event_rows = cursor.fetchall()
        for row in event_rows:
            insertion_time.append(row["insertion_time"])
            begin_time.append(row["begin_time"])
            reference_time.append(row["reference_time"])
            end_time.append(row["end_time"])
            # As each LYTAF files use same numbering for its event
            # types, change these numbers to distinguish different event
            # types from different LYTAF files.
            # lyra types -> 100s, e.g. lyra type 1 -> 101
            # manual types -> 200s, e.g. manual type 1 -> 201
            # ppt types -> 300s, e.g. ppt type 1 -> 301
            # science types -> 400s, e.g. science type 1 -> 401
            if suffix == "lyra":
                eventType_id.append(row["eventType_id"]+100)
            elif suffix == "ppt":
                eventType_id.append(row["eventType_id"]+300)
            elif suffix == "manual":
                eventType_id.append(row["eventType_id"]+200)
            elif suffix == "science":
                eventType_id.append(row["eventType_id"]+400)
        # Select and extract the event types from eventType table
        cursor.execute("SELECT * FROM eventType")
        eventType_rows = cursor.fetchall()
        for row in eventType_rows:
            _type.append(row["type"])
            _definition.append(row["definition"])
            # Alter id numbers so event types from different LYTAF
            # files can be distinguished, as above.
            if suffix == "lyra":
                _id.append(row["id"]+100)
            elif suffix == "ppt":
                _id.append(row["id"]+300)
            elif suffix == "manual":
                _id.append(row["id"]+200)
            elif suffix == "science":
                _id.append(row["id"]+400)
        
        # Sort arrays in order of increasing start time.
        # INSERT CODE HERE

        # Convert times to datetime objects
        # INSERT CODE HERE

        # Create dictionary of results and return it
        lytaf = {"insertion_time": insertion_time,
                 "begin_time": begin_time
                 "reference_time": reference_time
                 "end_time": end_time
                 "eventType_id": eventType_id
                 "type": np.asanyarray(_type, dtype=str)
                 "definition": np.asanyarray(_definition, dtype=str)
                 }

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
