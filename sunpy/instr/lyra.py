from __future__ import absolute_import, division, print_function

import os.path
import datetime
from warnings import warn
import copy
import csv
import sqlite3

import numpy as np
from astropy.io import fits
import pandas

from sunpy.time import parse_time
from sunpy import config
from sunpy.util.net import check_download_file
from sunpy import lightcurve

from sunpy.extern.six.moves import urllib

LYTAF_REMOTE_PATH = "http://proba2.oma.be/lyra/data/lytaf/"
LYTAF_PATH = config.get("downloads", "download_dir")


def remove_lytaf_events_from_lightcurve(lc, artifacts=None,
                                        return_artifacts=False,
                                        lytaf_path=None,
                                        force_use_local_lytaf=False):
    """
    Removes periods of LYRA artifacts defined in LYTAF from a LYRALightCurve.

    Parameters
    ----------
    lc : `sunpy.lightcurve.LightCurve`

     artifacts : list of strings
        Contain the artifact types to be removed.  For list of artifact types
        see reference [1].  For example, if user wants to remove only large
        angle rotations, listed at reference [1] as LAR, let artifacts=["LAR"].
        Default=[], i.e. no artifacts will be removed.

    return_artifacts : `bool`
        Set to True to return a `numpy.recarray` containing the start time, end
        time and type of all artifacts removed.
        Default=False

    lytaf_path : `str`
        directory path where the LYRA annotation files are stored.

    force_use_local_lytaf : `bool`
        Ensures current local version of lytaf files are not replaced by
        up-to-date online versions even if current local lytaf files do not
        cover entire input time range etc.
        Default=False

    Returns
    -------
    lc_new : `sunpy.lightcurve.LightCurve`
        copy of input LYRALightCurve with periods corresponding to artifacts
        removed.

    artifact_status : `dict`
        List of 4 variables containing information on what artifacts were
        found, removed, etc. from the time series.
        artifact_status["lytaf"] = artifacts found : `numpy.recarray`
            The full LYRA annotation file for the time series time range
            output by get_lytaf_events().
        artifact_status["removed"] = artifacts removed : `numpy.recarray`
            Artifacts which were found and removed from from time series.
        artifact_status["not_removed"] = artifacts found but not removed :
              `numpy.recarray`
            Artifacts which were found but not removed as they were not
            included when user defined artifacts kwarg.
        artifact_status["not_found"] = artifacts not found : `list` of strings
            Artifacts listed to be removed by user when defining
            artifacts kwarg which were not found in time series time range.

    References
    ----------
    [1] http://proba2.oma.be/data/TARDIS

    Examples
    --------
    Remove LARs (Large Angle Rotations) from LYRALightCurve for 4-Dec-2014:

        >>> import sunpy.lightcurve as lc
        >>> lc = lc.LYRALightCurve.create("2014-12-02")
        >>> lc_nolars = lc.remove_artifacts_from_lyralightcurve(lc, artifacts=["LAR"])

    To also retrieve information on the artifacts during that day:
        >>> lc_nolars, artifact_status = lc.remove_artifacts_from_lyralightcurve(
                lc, artifacts=["LAR"], return_artifacts=True)

    """
    # Check that input argument is of correct type
    if not lytaf_path:
        lytaf_path = LYTAF_PATH
    if not isinstance(lc, lightcurve.LightCurve):
        raise TypeError("lc must be a LightCurve object.")
    # Remove artifacts from time series
    data_columns = lc.data.columns
    time, channels, artifact_status = _remove_lytaf_events(
        lc.data.index,
        channels=[np.asanyarray(lc.data[col]) for col in data_columns],
        artifacts=artifacts, return_artifacts=True, lytaf_path=lytaf_path,
        force_use_local_lytaf=force_use_local_lytaf)
    # Create new copy copy of lightcurve and replace data with
    # artifact-free time series.
    lc_new = copy.deepcopy(lc)
    lc_new.data = pandas.DataFrame(
        index=time, data=dict((col, channels[i])
                              for i, col in enumerate(data_columns)))
    if return_artifacts:
        return lc_new, artifact_status
    else:
        return lc_new

def _remove_lytaf_events(time, channels=None, artifacts=None,
                         return_artifacts=False, fitsfile=None,
                         csvfile=None, filecolumns=None,
                         lytaf_path=None, force_use_local_lytaf=False):
    """
    Removes periods of LYRA artifacts from a time series.

    This functions removes periods corresponding to certain artifacts recorded
    in the LYRA annotation file from an array of times given by the time input.
    If a list of arrays of other properties is supplied through the channels
    kwarg, then the relevant values from these arrays are also removed.  This
    is done by assuming that each element in each array supplied corresponds to
    the time in the same index in time array.  The artifacts to be removed are
    given via the artifacts kwarg.  The default is "all", meaning that all
    artifacts will be removed.  However, a subset of artifacts can be removed
    by supplying a list of strings of the desired artifact types.

    Parameters
    ----------
    time : `numpy.ndarray` of `datetime.datetime`
        Gives the times of the timeseries.

    channels : `list` of `numpy.array` convertible to float64.
        Contains arrays of the irradiances taken at the times in the time
        variable.  Each element in the list must have the same number of
        elements as time.

    artifacts : `list` of strings
        Contain the artifact types to be removed.  For list of artifact types
        see reference [1].  For example, if user wants to remove only large
        angle rotations, listed at reference [1] as LAR, let artifacts=["LAR"].
        Default=[], i.e. no artifacts will be removed.

    return_artifacts : `bool`
        Set to True to return a numpy recarray containing the start time, end
        time and type of all artifacts removed.
        Default=False

    fitsfile : `str`
        file name (including file path and suffix, .fits) of output fits file
        which is generated if this kwarg is not None.
        Default=None, i.e. no fits file is output.

    csvfile : `str`
        file name (including file path and suffix, .csv) of output csv file
        which is generated if this kwarg is not None.
        Default=None, i.e. no csv file is output.

    filecolumns : `list` of strings
        Gives names of columns of any output files produced.  Although
        initially set to None above, the default is in fact
        ["time", "channel0", "channel1",..."channelN"]
        where N is the number of irradiance arrays in the channels input
        (assuming 0-indexed counting).

    lytaf_path : `str`
        directory path where the LYRA annotation files are stored.

    force_use_local_lytaf : `bool`
        Ensures current local version of lytaf files are not replaced by
        up-to-date online versions even if current local lytaf files do not
        cover entire input time range etc.
        Default=False

    Returns
    -------
    clean_time : `numpy.ndarray` of `datetime.datetime`
        time array with artifact periods removed.

    clean_channels : `list` ndarrays/array-likes convertible to float64
        list of irradiance arrays with artifact periods removed.

    artifact_status : `dict`
        List of 4 variables containing information on what artifacts were
        found, removed, etc. from the time series.
        artifact_status["lytaf"] = artifacts found : `numpy.recarray`
            The full LYRA annotation file for the time series time range
            output by get_lytaf_events().
        artifact_status["removed"] = artifacts removed : `numpy.recarray`
            Artifacts which were found and removed from from time series.
        artifact_status["not_removed"] = artifacts found but not removed :
              `numpy.recarray`
            Artifacts which were found but not removed as they were not
            included when user defined artifacts kwarg.
        artifact_status["not_found"] = artifacts not found : `list` of strings
            Artifacts listed to be removed by user when defining artifacts
            kwarg which were not found in time series time range.

    References
    ----------
    [1] http://proba2.oma.be/data/TARDIS

    Example
    -------
    Sample data for example
        >>> time = np.array([datetime(2013, 2, 1)+timedelta(minutes=i)
                             for i in range(120)])
        >>> channel_1 = np.zeros(len(TIME))+0.4
        >>> channel_2 = np.zeros(len(TIME))+0.1
    Remove LARs (Large Angle Rotations) from time series.
        >>> time_clean, channels_clean = remove_lyra_artifacts(
              time, channels=[channel_1, channel2], artifacts=['LAR'])

    """
    # Check inputs
    if not lytaf_path:
        lytaf_path = LYTAF_PATH
    if channels and type(channels) is not list:
        raise TypeError("channels must be None or a list of numpy arrays "
                        "of dtype 'float64'.")
    if not artifacts:
        raise ValueError("User has supplied no artifacts to remove.")
    if type(artifacts) is str:
        artifacts = [artifacts]
    if not all(isinstance(artifact_type, str) for artifact_type in artifacts):
        raise TypeError("All elements in artifacts must in strings.")
    all_lytaf_event_types = get_lytaf_event_types(lytaf_path=lytaf_path,
                                                  print_event_types=False)
    for artifact in artifacts:
        if artifact not in all_lytaf_event_types:
            print(all_lytaf_event_types)
            raise ValueError("{0} is not a valid artifact type. See above.".format(artifact))
    # Define outputs
    clean_time = np.array([parse_time(t) for t in time])
    clean_channels = copy.deepcopy(channels)
    artifacts_not_found = []
    # Get LYTAF file for given time range
    lytaf = get_lytaf_events(time[0], time[-1], lytaf_path=lytaf_path,
                             force_use_local_lytaf=force_use_local_lytaf)

    # Find events in lytaf which are to be removed from time series.
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
    if not len(artifact_indices):
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
        clean_time = np.delete(clean_time, bad_indices)
        if channels:
            for i, f in enumerate(clean_channels):
                clean_channels[i] = np.delete(f, bad_indices)
    # If return_artifacts kwarg is True, return a list containing
    # information on what artifacts found, removed, etc.  See docstring.
    if return_artifacts:
        artifact_status = {"lytaf": lytaf,
                           "removed": lytaf[artifact_indices],
                           "not_removed": np.delete(lytaf, artifact_indices),
                           "not_found": artifacts_not_found}
    # Output FITS file if fits kwarg is set
    if fitsfile:
        # Create time array of time strings rather than datetime objects
        # and verify filecolumns have been correctly input.  If None,
        # generate generic filecolumns (see docstring of function called
        # below.
        string_time, filecolumns = _prep_columns(time, channels, filecolumns)
        # Prepare column objects.
        cols = [fits.Column(name=filecolumns[0], format="26A",
                            array=string_time)]
        if channels:
            for i, f in enumerate(channels):
                cols.append(fits.Column(name=filecolumns[i+1], format="D",
                                        array=f))
        coldefs = fits.ColDefs(cols)
        tbhdu = fits.new_table(coldefs)
        hdu = fits.PrimaryHDU()
        tbhdulist = fits.HDUList([hdu, tbhdu])
        # Write data to fits file.
        tbhdulist.writeto(fitsfile)
    # Output csv file if csv kwarg is set.
    if csvfile:
        # Create time array of time strings rather than datetime objects
        # and verify filecolumns have been correctly input.  If None,
        # generate generic filecolumns (see docstring of function called
        # below.
        string_time, filecolumns = _prep_columns(time, channels, filecolumns)
        # Open and write data to csv file.
        with open(csvfile, 'w') as openfile:
            csvwriter = csv.writer(openfile, delimiter=';')
            # Write header.
            csvwriter.writerow(filecolumns)
            # Write data.
            if not channels:
                for i in range(len(time)):
                    csvwriter.writerow(string_time[i])
            else:
                for i in range(len(time)):
                    row = [string_time[i]]
                    for f in channels:
                        row.append(f[i])
                    csvwriter.writerow(row)
    # Return values.
    if return_artifacts:
        if not channels:
            return clean_time, artifact_status
        else:
            return clean_time, clean_channels, artifact_status
    else:
        if not channels:
            return clean_time
        else:
            return clean_time, clean_channels


def get_lytaf_events(start_time, end_time, lytaf_path=None,
                     combine_files=("lyra", "manual", "ppt", "science"),
                     csvfile=None, force_use_local_lytaf=False):
    """
    Extracts combined lytaf file for given time range.

    Given a time range defined by start_time and end_time, this function
    extracts the segments of each LYRA annotation file and combines them.

    Parameters
    ----------
    start_time : `datetime.datetime` or `str`
        Start time of period for which annotation file is required.

    end_time : `datetime.datetime` or `str`
        End time of period for which annotation file is required.

    lytaf_path : `str`
        directory path where the LYRA annotation files are stored.

    combine_files : `tuple` of strings
        States which LYRA annotation files are to be combined.
        Default is all four, i.e. lyra, manual, ppt, science.
        See Notes section for an explanation of each.

    force_use_local_lytaf : `bool`
        Ensures current local version of lytaf files are not replaced by
        up-to-date online versions even if current local lytaf files do not
        cover entire input time range etc.
        Default=False

    Returns
    -------
    lytaf : `numpy.recarray`
        Containing the various parameters stored in the LYTAF files.

    Notes
    -----
    There are four LYRA annotation files which mark different types of events
    or artifacts in the data.  They are named annotation_suffix.db where
    suffix is a variable equalling either lyra, manual, ppt, or science.

    annotation_lyra.db : contains entries regarding possible effects to
        the data due to normal operation of LYRA instrument.

    annotation_manual.db : contains entries regarding possible effects
        to the data due to unusual or manually logged events.

    annotation_ppt.db : contains entries regarding possible effects to
        the data due to pointing or positioning of PROBA2.

    annotation_science.db : contains events in the data scientifically
        interesting, e.g. GOES flares.

    References
    ----------
    Further documentation: http://proba2.oma.be/data/TARDIS

    Examples
    --------
    Get all events in the LYTAF files for January 2014
        >>> from sunpy.instr.lyra import get_lytaf_events
        >>> lytaf = get_lytaf_events('2014-01-01', '2014-02-01')

    """
    # Check inputs
    # Check lytaf path
    if not lytaf_path:
        lytaf_path = LYTAF_PATH
    # Check start_time and end_time is a date string or datetime object
    start_time = parse_time(start_time)
    end_time = parse_time(end_time)
    # Check combine_files contains correct inputs
    if not all(suffix in ["lyra", "manual", "ppt", "science"]
               for suffix in combine_files):
        raise ValueError("Elements in combine_files must be strings equalling "
                         "'lyra', 'manual', 'ppt', or 'science'.")
    # Remove any duplicates from combine_files input
    combine_files = list(set(combine_files))
    combine_files.sort()
    # Convert input times to UNIX timestamp format since this is the
    # time format in the annotation files
    start_time_uts = (start_time - datetime.datetime(1970, 1, 1)).total_seconds()
    end_time_uts = (end_time - datetime.datetime(1970, 1, 1)).total_seconds()

    # Define numpy record array which will hold the information from
    # the annotation file.
    lytaf = np.empty((0,), dtype=[("insertion_time", object),
                                  ("begin_time", object),
                                  ("reference_time", object),
                                  ("end_time", object),
                                  ("event_type", object),
                                  ("event_definition", object)])
    # Access annotation files
    for suffix in combine_files:
        # Check database files are present
        dbname = "annotation_{0}.db".format(suffix)
        check_download_file(dbname, LYTAF_REMOTE_PATH, lytaf_path)
        # Open SQLITE3 annotation files
        connection = sqlite3.connect(os.path.join(lytaf_path, dbname))
        # Create cursor to manipulate data in annotation file
        cursor = connection.cursor()
        # Check if lytaf file spans the start and end times defined by
        # user.  If not, download newest version.
        # First get start time of first event and end time of last
        # event in lytaf.
        cursor.execute("select begin_time from event order by begin_time asc "
                       "limit 1;")
        db_first_begin_time = cursor.fetchone()[0]
        db_first_begin_time = datetime.datetime.fromtimestamp(db_first_begin_time)
        cursor.execute("select end_time from event order by end_time desc "
                       "limit 1;")
        db_last_end_time = cursor.fetchone()[0]
        db_last_end_time = datetime.datetime.fromtimestamp(db_last_end_time)
        # If lytaf does not include entire input time range...
        if not force_use_local_lytaf:
            if end_time > db_last_end_time or start_time < db_first_begin_time:
                # ...close lytaf file...
                cursor.close()
                connection.close()
                # ...Download latest lytaf file...
                check_download_file(dbname, LYTAF_REMOTE_PATH, lytaf_path,
                                    replace=True)
                # ...and open new version of lytaf database.
                connection = sqlite3.connect(os.path.join(lytaf_path, dbname))
                cursor = connection.cursor()
        # Select and extract the data from event table within file within
        # given time range
        cursor.execute("select insertion_time, begin_time, reference_time, "
                       "end_time, eventType_id from event where end_time >= "
                       "{0} and begin_time <= "
                       "{1}".format(start_time_uts, end_time_uts))
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
                              np.array((datetime.datetime.utcfromtimestamp(event_row[0]),
                                        datetime.datetime.utcfromtimestamp(event_row[1]),
                                        datetime.datetime.utcfromtimestamp(event_row[2]),
                                        datetime.datetime.utcfromtimestamp(event_row[3]),
                                        eventType_type[id_index],
                                        eventType_definition[id_index]), dtype=lytaf.dtype))
        # Close file
        cursor.close()
        connection.close()
    # Sort lytaf in ascending order of begin time
    np.recarray.sort(lytaf, order="begin_time")

    # If csvfile kwarg is set, write out lytaf to csv file
    if csvfile:
        # Open and write data to csv file.
        with open(csvfile, 'w') as openfile:
            csvwriter = csv.writer(openfile, delimiter=';')
            # Write header.
            csvwriter.writerow(lytaf.dtype.names)
            # Write data.
            for row in lytaf:
                new_row = []
                new_row.append(row[0].strftime("%Y-%m-%dT%H:%M:%S"))
                new_row.append(row[1].strftime("%Y-%m-%dT%H:%M:%S"))
                new_row.append(row[2].strftime("%Y-%m-%dT%H:%M:%S"))
                new_row.append(row[3].strftime("%Y-%m-%dT%H:%M:%S"))
                new_row.append(row[4])
                new_row.append(row[5])
                csvwriter.writerow(new_row)

    return lytaf

def get_lytaf_event_types(lytaf_path=None, print_event_types=True):
    """Prints the different event types in the each of the LYTAF databases.

    Parameters
    ----------
    lytaf_path : `str`
        Path location where LYTAF files are stored.
        Default = LYTAF_PATH defined above.

    print_event_types : `bool`
        If True, prints the artifacts in each lytaf database to screen.

    Returns
    -------
    all_event_types : `list`
        List of all events types in all lytaf databases.

    """
    # Set lytaf_path is not done by user
    if not lytaf_path:
        lytaf_path = LYTAF_PATH
    suffixes = ["lyra", "manual", "ppt", "science"]
    all_event_types = []
    # For each database file extract the event types and print them.
    if print_event_types:
        print("\nLYTAF Event Types\n-----------------\n")
    for suffix in suffixes:
        dbname = "annotation_{0}.db".format(suffix)
        # Check database file exists, else download it.
        check_download_file(dbname, LYTAF_REMOTE_PATH, lytaf_path)
        # Open SQLITE3 LYTAF files
        connection = sqlite3.connect(os.path.join(lytaf_path, dbname))
        # Create cursor to manipulate data in annotation file
        cursor = connection.cursor()
        cursor.execute("select type from eventType;")
        event_types = cursor.fetchall()
        all_event_types.append(event_types)
        if print_event_types:
            print("----------------\n{0} database\n----------------"
                  .format(suffix))
            for event_type in event_types:
                print(str(event_type[0]))
            print(" ")
    # Unpack event types in all_event_types into single list
    all_event_types = [event_type[0] for event_types in all_event_types
                       for event_type in event_types]
    return all_event_types


def download_lytaf_database(lytaf_dir=''):
    """download latest Proba2 pointing database from Proba2 Science Center"""
    url = 'http://proba2.oma.be/lyra/data/lytaf/annotation_ppt.db'
    destination = os.path.join(lytaf_dir, 'annotation_ppt.db')
    urllib.request.urlretrieve(url, destination)

    return


def split_series_using_lytaf(timearray, data, lytaf):
    """
    Proba-2 analysis code for splitting up LYRA timeseries around locations
    where LARs (and other data events) are observed.

    Parameters
    ----------
    timearray : `numpy.ndarray` of times understood by `sunpy.time.parse_time`
        function.
    data : `numpy.array` corresponding to the given time array
    lytaf : `numpy.recarray`
        Events obtained from querying LYTAF database using
        lyra.get_lytaf_events().

    Output
    ------
    output : `list` of dictionaries
        Each dictionary contains a sub-series corresponding to an interval of
        'good data'.
    """
    n = len(timearray)
    mask = np.ones(n)
    el = len(lytaf)

    # make the input time array a list of datetime objects
    datetime_array = []
    for tim in timearray:
        datetime_array.append(parse_time(tim))

    # scan through each entry retrieved from the LYTAF database
    for j in range(0, el):
        # want to mark all times with events as bad in the mask, i.e. = 0
        start_dt = lytaf['begin_time'][j]
        end_dt = lytaf['end_time'][j]

        # find the start and end indices for each event
        start_ind = np.searchsorted(datetime_array, start_dt)
        end_ind = np.searchsorted(datetime_array, end_dt)

        # append the mask to mark event as 'bad'
        mask[start_ind:end_ind] = 0

    diffmask = np.diff(mask)
    tmp_discontinuity = np.where(diffmask != 0.)
    # disc contains the indices of mask where there are discontinuities
    disc = tmp_discontinuity[0]

    if len(disc) == 0:
        print('No events found within time series interval. '
              'Returning original series.')
        return [{'subtimes': datetime_array, 'subdata': data}]

    # -1 in diffmask means went from good data to bad
    # +1 means went from bad data to good

    # want to get the data between a +1 and the next -1

    # if the first discontinuity is a -1 then the start of the series was good.
    if diffmask[disc[0]] == -1.0:
        # make sure we can always start from disc[0] below
        disc = np.insert(disc, 0, 0)

    split_series = []

    limit = len(disc)
    # now extract the good data regions and ignore the bad ones
    for h in range(0, limit, 2):

        if h == limit-1:
            # can't index h+1 here. Go to end of series
            subtimes = datetime_array[disc[h]:-1]
            subdata = data[disc[h]:-1]
            subseries = {'subtimes':subtimes, 'subdata':subdata}
            split_series.append(subseries)
        else:
            subtimes = datetime_array[disc[h]:disc[h+1]]
            subdata = data[disc[h]:disc[h+1]]
            subseries = {'subtimes':subtimes, 'subdata':subdata}
            split_series.append(subseries)

    return split_series

def _lytaf_event2string(integers):
    if type(integers) == int:
        integers = [integers]
    #else:
    #    n=len(integers)
    out = []

    for i in integers:
        if i == 1:
            out.append('LAR')
        if i == 2:
            out.append('N/A')
        if i == 3:
            out.append('UV occult.')
        if i == 4:
            out.append('Vis. occult.')
        if i == 5:
            out.append('Offpoint')
        if i == 6:
            out.append('SAA')
        if i == 7:
            out.append('Auroral zone')
        if i == 8:
            out.append('Moon in LYRA')
        if i == 9:
            out.append('Moon in SWAP')
        if i == 10:
            out.append('Venus in LYRA')
        if i == 11:
            out.append('Venus in SWAP')

    return out

def _prep_columns(time, channels=None, filecolumns=None):
    """
    Checks and prepares data to be written out to a file.

    Firstly, this function converts the elements of time, whose entries are
    assumed to be datetime objects, to time strings.  Secondly, it checks
    whether the number of elements in an input list of column names,
    filecolumns, is equal to the number of arrays in the list, channels.
    If not, a ValueError is raised.  If however filecolumns equals None, a
    filenames list is generated equal to ["time", "channel0", "channel1",...,
    "channelN"] where N is the number of arrays in the list, channels
    (assuming 0-indexed counting).

    """
    # Convert time which contains datetime objects to time strings.
    string_time = np.array([t.strftime("%Y-%m-%dT%H:%M:%S.%f") for t in time])
    # If filenames is given...
    if filecolumns:
        # ...check all the elements are strings...
        if all(isinstance(column, str) for column in filecolumns) is False:
            raise TypeError("All elements in filecolumns must by strings.")
        # ...and that there are the same number of elements as there
        # are arrays in channels, plus 1 for a time array.  Otherwise
        # raise a ValueError.
        if channels:
            ncol = 1 + len(channels)
        else:
            ncol = 1
        if len(filecolumns) != ncol:
            raise ValueError("Number of elements in filecolumns must be "
                             "equal to the number of input data arrays, "
                             "i.e. time + elements in channels.")
    # If filenames not given, create a list of columns names of the
    # form: ["time", "channel0", "channel1",...,"channelN"] where N
    # is the number of arrays in channels (assuming 0-indexed counting).
    else:
        if channels:
            filecolumns = ["channel{0}".format(fluxnum)
                           for fluxnum in range(len(channels))]
            filecolumns.insert(0, "time")
        else:
            filecolumns = ["time"]

    return string_time, filecolumns
