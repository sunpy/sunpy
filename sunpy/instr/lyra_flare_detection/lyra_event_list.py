from __future__ import absolute_import
from __future__ import division

import os.path
from datetime import datetime
from warnings import warn
import copy
import csv
import urllib

import numpy as np
import sqlite3
from itertools import chain
from sunpy.time import parse_time
from astropy.io import fits

from datetime import timedelta
import matplotlib.pyplot as plt

RISE_FACTOR = 1.005
FALL_FACTOR = 0.5
# Set mean daily minimum irradiance in Zr channel from first light
# (Jan 2010) until mid 2014.
NORM = 0.001

def lyra_event_list(tstart, tend,
                    data_path=os.path.expanduser(
                        os.path.join("~", "pro", "data", "LYRA", "fits"))):
    """
    Returns a list of LYRA flares based on input start and end time.

    """

def get_lyra_timeseries(tstart, tend, level=2, unit_preference=["std", "bst"],
                        data_path=os.path.expanduser(
                            os.path.join("~", "pro", "data", "LYRA", "fits"))):
    """
    Returns LYRA timeseries for a given time range.

    This function returns a numpy record array containing the timeseries of
    the four LYRA channels between the input start and end time.  This is

    Parameters
    ----------
    tstart : datetime object or timestring
        Start time of time range.

    tend : datetime object or timestring
        End time of time range.

    level : int
        Level of LYRA data user wishes to use.  Options are 1, 2, or 3.
        See reference [1] for explanations on the differences between each
        level of data.

    unit_preference : list of strings
        Denotes whether user prefers data from nominal or back-up units.
        'std' denotes nominal unit, 'bst' denotes back-up units.
        ['std', 'bst'] means data from nominal unit will be searched for
        first, but if such data can't be found, back-up unit data will be
        used.
        ['bst', 'std'] means data from back-up units will be searched for
        first, but if such data can't be found, nominal unit data will be
        used.
        ['std'] means only data from nominal unit will be searched for.  If
        nominal unit data cannot be found, a data gap will appear in
        timeseries.
        ['bst'] means only data from back-up unit will be searched for.  If
        back-up unit data cannot be found, a data gap will appear in
        timeseries.

    data_path : string
        directory path in which to save/search for required LYRA fits files.

    Returns
    -------
    timeseries : numpy record array
        Resulting timeseries for the time range in question.

    References
    ----------
    [1]  http://proba2.oma.be/data/LYRA

    Examples
    --------

    """
    # Check inputs are of correct types etc.
    parse_time(tstart)
    if type(tstart) is not datetime.datetime:
        raise TypeError("tstart must be a datetime object or valid time "
                        "string.")
    parse_time(tend)
    if type(tend) is not datetime.datetime:
        raise TypeError("tend must be a datetime object or valid time string.")
    if tstart >= tend:
        raise ValueError("Start time must be before end time.")
    if type(level) is not int:
        raise TypeError("{0} must be an int.".format(level))
    elif level != 1 or level != 2 or level != 3:
        raise ValueError("{0} must be equal to either 1, 2, or "
                         "3.".format(level))
    if type(unit_preference) != list:
        raise TypeError("unit_preference must be a list of strings.")
    if not all(unit=='std' or unit=='bst' for unit in unit_preference):
        raise ValueError("All elements in unit_perference must be either "
                         "'std' or 'bst'.")
    if len(unit_preference) < 1 or len(unit_preference) > 2:
        raise ValueError("unit_preference must have either 1 or 2 elements.")
    if not os.path.isdir(data_path):
        raise ValueError("{0} is not a valid directory.".format(data_path))
    # Define values needed later
    num_units = len(unit_preference)
    remote_level_dir = ['eng', 'bsd', 'bsd']
    
    # Determine what files are needed for requested time range.
    # Find dates of start and end times.
    filedates = [datetime(tstart.year, tstart.month, tstart.day),
                 datetime(tend.year, tend.month, tend.day)]
    # Create list of LYRA fits file name tags for each required date
    # and enter date of start time.
    filetags = ["lyra_{0}{1}{2}-000000_lev{3}".format(
        filedates[0].year, filedates[0].month, filedates[0].day, level)]
    if filedates[0] != filedates[-1]:
        # If start and end time are on different days create tags for
        # days after start end time
        dayrange = range(1, (filedates[1]-filedates[0]).days)
        for i in dayrange:
            filedate = filedates[0]+timedelta(i)
            filetags.append("lyra_{0}{1}{2}-000000_lev{3}".format(
                filedate.year, filedate.month, filedate.day, level))

    # Search if files exist in data_path.  If not, download them.
    filenames = []
    for tag in filetags:
        filenames.append(os.path.join(
            data_path, tag, "_{0}.fits".format(unit_preference[0])))
        if not os.path.isfile(filenames[-1]):
            getfile = urllib.URLopener()
            try:
                getfile.retrieve(
                    "http://http://proba2.oma.be/lyra/data/{0}/{1}/{2}/{3}/"
                    "{4}_{5}.fits".format(
                        remote_level_dir[level-1], tag[5:9], tag[9:11],
                        tag[11:13], tag, unit_preference[0]), filenames[-1])
            except IOError:
                # If data file for first preference unit cannot be found
                # Try downloading data file for second preference unit.
                if num_units > 1:
                    try:
                        getfile.retrieve(
                            "http://http://proba2.oma.be/lyra/data/{0}/{1}/"
                            "{2}/{3}/{4}_{5}.fits".format(
                                remote_level_dir[level-1], tag[5:9], tag[9:11],
                                tag[11:13], tag, unit_preference[0]),
                            filenames[-1])
                    except IOError:
                        del(filenames[-1])
                        raise Warning("File for {0}/{1}/{2} could not be "
                                      "found.".format(tag[5:9], tag[9:11],
                                                      tag[11:13]))
                else:
                    del(filenames[-1])
                    raise Warning("File for {0}/{1}/{2} could not be "
                                      "found.".format(tag[5:9], tag[9:11],
                                                      tag[11:13]))
    # 
            
    
    

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
    time = _check_datetime(time)
    # Define variables to be used later
    flare_indices = []
    lyra_events = np.empty((0,), dtype=[("start_time", object),
                                        ("peak_time", object),
                                        ("end_time", object),
                                        ("comments", object)])
    # object LYRA artifacts from timeseries
    clean_time, fluxlist, artifact_status = remove_lyra_artifacts(time, [flux],
        artifacts=["UV occ.", "Offpoint", "LAR", "Calibration", "SAA",
                   "Vis occ.", "Operational Anomaly", "Glitch", "ASIC reload"],
                   return_artifacts=True)
    clean_flux = fluxlist[0]
    artifacts_removed = artifact_status[1]
    # Perform subtraction so median irradiance of time series is at
    # average daily minimum from first 4 years of mission.
    clean_flux = clean_flux - (np.median(clean_flux)-NORM)
    # Get derivative of flux wrt time
    time_timedelta = clean_time[1:-1]-clean_time[0:-2]
    dt = np.zeros(len(time_timedelta), dtype="float64")
    for i, t, in enumerate(time_timedelta):
        dt[i] = t.total_seconds()
    dfdt = np.gradient(clean_flux[0:-2], dt)
    # Get locations where derivative is positive
    pos_deriv = np.where(dfdt > 0)[0]
    neg_deriv = np.where(dfdt < 0)[0]
    pos_deriv0 = np.where(dfdt >= 0)[0]
    neg_deriv0 = np.where(dfdt <= 0)[0]
    # Find difference between each time point and the one 4
    # observations ahead.
    time_timedelta4 = clean_time[4:-1]-clean_time[0:-5]
    dt4 = np.zeros(len(time_timedelta4))
    for i, t, in enumerate(time_timedelta4):
        dt4[i] = t.total_seconds()
    # Find all possible flare start times.
    i=0
    while i < len(pos_deriv)-4:
        # Start time criteria
        if (pos_deriv[i:i+4]-pos_deriv[i] == np.arange(4)).all() and \
          dt4[pos_deriv[i]] > 210 and dt4[pos_deriv[i]] < 270 and \
          clean_flux[pos_deriv[i+4]]/clean_flux[pos_deriv[i]] >= RISE_FACTOR:
            # Find start time which is defined as earliest continuous
            # increase in flux before the point found by the above
            # criteria.
            try:
                k = np.where(neg_deriv0 < pos_deriv[i])[0][-1]
                kk = np.where(pos_deriv > neg_deriv0[k])[0][0]
            except IndexError:
                kk = i                
            start_index = pos_deriv[kk]
            # If artifact is at start of flare, set start time to
            # directly afterwards.
            artifact_check = np.logical_and(
                artifacts_removed["end_time"] > clean_time[start_index],
                artifacts_removed["end_time"] < clean_time[pos_deriv[kk+2]])
            if artifact_check.any() == True:
                new_index = \
                  np.where(clean_time[pos_deriv] > \
                           artifacts_removed[artifact_check][-1]["end_time"])
                start_index = pos_deriv[new_index[0][0]]
            # Next, find index of flare end time.
            jj = np.where(neg_deriv > start_index)[0][0]
            j = neg_deriv[jj]
            try:
                l = np.where(
                    clean_flux[j:] <= max(clean_flux[pos_deriv[i]:j]) - \
                    (max(clean_flux[pos_deriv[i]:j]) - \
                     clean_flux[pos_deriv[i]]) * FALL_FACTOR)[0][0]+j
            except IndexError:
                i = i+1
            else:
                try:
                    m = np.where(pos_deriv0 > l)[0][0]
                except IndexError:
                    i = i+1
                else:
                    end_index = pos_deriv0[m]-1
                    # If artifact is at end of flare, set end time to
                    # directly beforehand.
                    artifact_check = np.logical_and(
                      artifacts_removed["begin_time"] < clean_time[end_index],
                        artifacts_removed["begin_time"] > \
                        clean_time[end_index-2])
                    if artifact_check.any() == True:
                        new_index = np.where(clean_time < \
                            artifacts_removed[artifact_check][0]["begin_time"])
                        end_index = new_index[0][-1]
                    # find index of peak time
                    peak_index = np.where(
                        clean_flux == max(clean_flux[start_index:end_index]))
                    peak_index = peak_index[0][0]
                    # Record flare start, peak and end times
                    flare_indices.append((start_index, peak_index, end_index))
                    lyra_events = np.append(
                        lyra_events, np.empty(1, dtype=lyra_events.dtype))
                    lyra_events[-1]["start_time"] = clean_time[start_index]
                    lyra_events[-1]["peak_time"] = clean_time[peak_index]
                    lyra_events[-1]["end_time"] = clean_time[end_index]
                    # If the most recently found flare is during the decay phase
                    # of another reset end time of previous flare to start time
                    # of this flare.
                    if len(lyra_events) > 1:
                        if lyra_events[-2]["end_time"] > \
                          lyra_events[-1]["start_time"]:
                            lyra_events[-2]["end_time"] = \
                              lyra_events[-1]["start_time"]
                    # Finally, set principle iterator, i, to the peak of the
                    # flare just found so that algorithm will start looking
                    # for flares during the decay phase of this flare and
                    # beyond.  This ensures that flares during the decay
                    # phase are also located.
                    i = np.where(pos_deriv > peak_index)[0][0]
        else:
            i = i+1

    return lyra_events

def remove_lyra_artifacts(time, fluxes=None, artifacts="All",
                          return_artifacts=False, fitsfile=None,
                          csvfile=None, filecolumns=None):
    """
    Removes periods of LYRA artifacts from a time series.

    This functions removes periods correspoding to certain artifacts recorded
    in the LYRA annotation file from an array of times given by the time input.
    If a list of arrays of other properties is supplied through the fluxes
    kwarg, then the relevant values from these arrays are also removed.  This
    is done by assuming that each element in each array supplied coresponds to
    the time in the same index in time array.  The artifacts to be removed are
    given via the artifacts kwarg.  The default is "all", meaning that all
    artifacts will be removed.  However, a subset of artifacts can be removed
    by supplying a list of strings of the desired artifact types.

    Parameters
    ----------
    time : ndarray/array-like of datetime objects
        Gives the times of the timeseries.

    fluxes : (optional) list of ndarrays/array-likes convertible to float64.
        Contains the fluxes/properties taken at the times in the time array.
        Each element in the list must have the same number of elements as time.

    artifacts : list of strings
        Contain the artifact types to be removed.  For list of artifact types
        see reference [1].  For example, if user wants to remove only large
        angle rotations, listed at reference [1] as LAR, let artifacts=["LAR"].

    return_artifacts : (optional) bool
        Set to True to return a numpy recarray containing the start time, end
        time and type of all artifacts removed.
        Default=False

    fitsfile : (optional) string
        file name (including file path) of output fits file which is generated
        if this kwarg is not None.
        Default=None, i.e. not csv file is output.

    csvfile : (optional) string
        file name (including file path) of output csv file which is generated
        if this kwarg is not None.
        Default=None, i.e. not csv file is output.

    filecolumns : (optional) list of strings
        Gives names of columns of any output files produced.  Although
        initially set to None above, the default is in fact
        ["time", "fluxes0", "fluxes1",..."fluxesN"]
        where N is the number of flux arrays in the fluxes input
        (assuming 0-indexed counting).
        
    Returns
    -------
    clean_time : ndarray/array-like of datetime objects
        time array with artifact periods removed.

    clean_fluxes : (optional) list ndarrays/array-likes convertible to float64
        list of fluxes with artifact periods removed.

    artifacts_removed : (optional) numpy recarray
        Three columns, "start_time", "end_time", "type", containing start
        time, end time and type of artifacts removed.

    References
    ----------
    [1] http://proba2.oma.be/data/TARDIS

    Example
    -------

    """
    # Check inputs
    time = _check_datetime(time)
    if not all(isinstance(artifact_type, str) for artifact_type in artifacts):
        raise TypeError("All elements in artifacts must in strings.")
    if type(fluxes) is not None and type(fluxes) is not list:
        raise TypeError("fluxes must be None or a list of numpy arrays of "
                        "dtype 'float64'.")
    # Define outputs
    clean_time = copy.deepcopy(time)
    clean_fluxes = copy.deepcopy(fluxes)
    artifacts_not_found =[]
    # Get LYTAF file for given time range
    lytaf = extract_combined_lytaf(time[0], time[-1])
    
    # Find events in lytaf which are to be removed from time series.
    if artifacts == "All":
        artifact_indices = np.arange(len(lytaf["begin_time"]))
    else:
        artifact_indices = np.empty(0, dtype="int64")
        for artifact_type in artifacts:
            indices = np.where(lytaf["event_type"] == artifact_type)[0]
            # If none of a given type of artifact is found, record this
            # type in artifact_not_found list.
            if len(indices) == 0:
                artifacts_not_found.append(artifact_type)
            else:
                # Else, record the indices of the artifacts of this type
                artifact_indices = np.concatenate((artifact_indices, indices))
        artifact_indices.sort()

    # Remove relevant artifacts from timeseries. If none of the
    # artifacts the user wanted removed were found, raise a warning and
    # continue with code.
    if len(artifact_indices) == 0:
        warn("None of user supplied artifacts were found.")
        artifacts_not_found = artifacts
    else:
        # Remove periods corresponding to artifacts from flux and time
        # arrays.
        bad_indices = np.empty(0, dtype="int64")
        all_indices = np.arange(len(time))
        for index in artifact_indices:
            bad_period = np.logical_and(time >= lytaf["begin_time"][index],
                                        time <= lytaf["end_time"][index])
            bad_indices = np.append(bad_indices, all_indices[bad_period])
        clean_time = np.delete(time, bad_indices)
        if fluxes is not None:
            for i, f in enumerate(clean_fluxes):
                clean_fluxes[i] = np.delete(f, bad_indices)
    # If return_artifacts kwarg is True, return a list containing
    # information on what artifacts found, removed, etc.  See docstring.
    if return_artifacts is True:
        # Define output list for artifact info
        artifact_status = []
        if artifacts_not_found == artifacts:
            # Artifacts found in annotation file
            artifact_status.append(lytaf)
            # Artifacts removed
            artifact_status.append(None)
            # Artifacts not removed
            artifact_status.append(None)
            # Artifacts not found
            artifact_status.append(artifacts_not_found)
        else:
            # Artifacts found in annotation file
            artifact_status.append(lytaf)
            # Artifacts removed            
            artifact_status.append(lytaf[artifact_indices])
            # Artifacts not removed
            artifact_status.append(np.delete(lytaf, artifact_indices)) 
            # Artifacts not found
            if artifacts == "All":
                artifact_status.append(None)
            else:
                artifact_status.append(artifacts_not_found)

    # Output FITS file if fits kwarg is set
    if fitsfile != None:
        # Create time array of time strings rather than datetime objects
        # and verify filecolumns have been correctly input.  If None,
        # generate generic filecolumns (see docstring og function called
        # below.
        string_time, filecolumns = _prep_columns(time, fluxes, filecolumns)
        # Prepare column objects.
        cols = [fits.Column(name=filecolumns[0], format="26A",
                            array=string_time)]
        if fluxes != None:
            for i, f in enumerate(fluxes):
                cols.append(fits.Column(name=filecolumns[i+1], format="D",
                                        array=f))
        coldefs = fits.ColDefs(cols)
        tbhdu = fits.new_table(coldefs)
        hdu = fits.PrimaryHDU()
        tbhdulist = fits.HDUList([hdu, tbhdu])
        # Write data to fits file.
        tbhdulist.writeto(fitsfile)
    # Output csv file if fits kwarg is set.
    if csvfile != None:
        # Create time array of time strings rather than datetime objects
        # and verify filecolumns have been correctly input.  If None,
        # generate generic filecolumns (see docstring og function called
        # below.
        string_time, filecolumns = prep_columns(time, fluxes, filecolumns)
        # Open and write data to csv file.
        with open(csvfile, 'w') as openfile:
            csvwriter = csv.writer(openfile, delimiter=';')
            # Write header.
            csvwriter.writerow(filecolumns)
            # Write data.
            if fluxes == None:
                for i in range(len(time)):
                    csvwriter.writerow(string_time[i])
            else:
                for i in range(len(time)):
                    row = [string_time[i]]
                    for f in fluxes:
                        row.append(f[i])
                    csvwriter.writerow(row)
    # Return values.
    if return_artifacts is True:
        if fluxes is None:
            return clean_time, artifact_status
        else:
            return clean_time, clean_fluxes, artifact_status
    else:
        if fluxes is None:
            return clean_time
        else:
            return clean_time, clean_fluxes

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
    tstart_uts = (tstart - datetime(1970, 1, 1)).total_seconds()
    tend_uts = (tend - datetime(1970, 1, 1)).total_seconds()

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
            lytaf = np.append(
                lytaf, np.array((datetime.utcfromtimestamp(event_row[0]),
                                 datetime.utcfromtimestamp(event_row[1]),
                                 datetime.utcfromtimestamp(event_row[2]),
                                 datetime.utcfromtimestamp(event_row[3]),
                                 eventType_type[id_index],
                                 eventType_definition[id_index]),
                                 dtype=lytaf.dtype))
        # Close file
        cursor.close()
        connection.close()
    # Sort lytaf in ascending order of begin time
    np.recarray.sort(lytaf, order="begin_time")    

    #return event_rows, eventType_rows
    return lytaf

def _check_datetime(time):
    """
    Checks or tries to convert input array to array of datetime objects.

    Returns input time array with elements as datetime objects or raises an
    TypeError if time not of valid format.  Input format can be anything
    convertible to datetime by datetime() function or any time string valid as
    an input to sunpy.time.parse_time().

    """
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
                raise TypeError("time must be an array or array-like of "
                                "datetime objects or valid time strings.")
        else:
            raise TypeError("time must be an array or array-like of "
                            "datetime objects or valid time strings.")
    return time

def _prep_columns(time, fluxes, filecolumns):
    """
    Checks and prepares data to be written out to a file.

    Firstly, this function converts the elements of time, whose entries are
    assumed to be datetime objects, to time strings.  Secondly, it checks
    whether the number of elements in an input list of columns names,
    filenames, is equal to the number of arrays in the list, fluxes.  If not,
    a Value Error is raised.  If however filecolumns equals None, a filenames
    list is generated equal to ["time", "fluxes0", "fluxes1",...,"fluxesN"]
    where N is the number of arrays in the list, fluxes
    (assuming 0-indexed counting).

    """
    # Convert time which contains datetime objects to time strings.
    string_time = np.empty(len(time), dtype="S26")
    for i, t in enumerate(time):
        string_time[i] = t.strftime("%Y-%m-%dT%H:%M:%S.%f")

    # If filenames is given...
    if filecolumns != None:
        # ...check all the elements are strings...
        if all(isinstance(column, str) for column in filecolumns) is False:
            raise TypeError("All elements in filecolumns must by strings.")
        # ...and that there are the same number of elements as there
        # are arrays in fluxes, plus 1 for a time array.  Otherwise
        # raise a ValueError.
        if fluxes != None:
            ncol = 1 + len(fluxes)
        else:
            ncol = 1
        if len(filecolumns) != ncol:
            raise ValueError("Number of elements in filecolumns must be "
                             "equal to the number of input data arrays, "
                             "i.e. time + elements in fluxes.")
    # If filenames not given, create a list of columns names of the
    # form: ["time", "fluxes0", "fluxes1",...,"fluxesN"] where N is the
    # number of arrays in fluxes (assuming 0-indexed counting).
    else:
        if fluxes != None:
            filecolumns = ["fluxes{0}".format(fluxnum)
                           for fluxnum in range(len(fluxes))]
            filecolumns.insert(0, "time")
        else:
            filecolumns = ["time"]

    return string_time, filecolumns

def testing_find_lyra_events(find_events=False):
    fitspath = "../data/LYRA/fits/"
    #fitsname = "lyra_20140621-000000_lev3_std.fits"
    #fitsname = "lyra_20140101-000000_lev3_std.fits"
    #fitsname = "lyra_20140102-000000_lev3_std.fits"
    #fitsname = "lyra_20140103-000000_lev3_std.fits"
    #fitsname = "lyra_20140601-000000_lev3_std.fits"
    #fitsname = "lyra_20140602-000000_lev3_std.fits"
    #fitsname = "lyra_20140603-000000_lev3_std.fits"
    #fitsname = "lyra_20140604-000000_lev3_std.fits"
    #fitsname = "lyra_20140605-000000_lev3_std.fits"
    #fitsname = "lyra_20140606-000000_lev3_std.fits"
    #fitsname = "lyra_20140607-000000_lev3_std.fits"
    #fitsname = "lyra_20100601-000000_lev3_std.fits"
    fitsname = "lyra_20100602-000000_lev3_std.fits"
    fitsname = "lyra_20100603-000000_lev3_std.fits"
    #fitsname = "lyra_20100604-000000_lev3_std.fits"
    #fitsname = "lyra_20111001-000000_lev3_std.fits"
    fitsfile = fitspath + fitsname
    ly = fits.open(fitsfile)
    t = ly[1].data["TIME"]
    orig_flux = ly[1].data["CHANNEL4"]
    orig_time = np.empty(len(t), dtype="object")
    for i, tt in enumerate(t):
        orig_time[i] = datetime(int(fitsname[5:9]), int(fitsname[9:11]),
                                int(fitsname[11:13]), 0, 0) + \
                                timedelta(0,0,0,0,int(tt))
    time = copy.deepcopy(orig_time)
    flux = copy.deepcopy(orig_flux)
    time, flux = remove_lyra_artifacts(time, [flux],
        artifacts=["UV occ.", "Offpoint", "LAR", "Calibration", "SAA",
                   "Vis occ.", "Operational Anomaly", "Glitch", "ASIC reload"])
    flux = flux[0]
    if find_events is False:
        return orig_time, orig_flux, time, flux
    else:
        ev = find_lyra_events(flux, time)
        ind = []
        for i in range(len(ev)):
            ind.append(np.where(time==ev[i]["start_time"])[0][0])
            ind.append(np.where(time==ev[i]["end_time"])[0][0])
        plt.ion()
        plt.plot(time, flux, 'o', color='blue')
        plt.plot(time[ind], flux[ind], 'o', color='red')
        return orig_time, orig_flux, time, flux, ev, ind
