# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
__authors__ = ["Alex Hamilton, Stuart Mumford"]
__email__ = "stuart@mumford.me.uk"

from sunpy.util.metadata import MetaDict
import itertools
import copy

import warnings
import inspect

from sunpy.time import TimeRange, parse_time


class TimeSeriesMetaData(object):
    """
    An object used to store metadata for TimeSeries objects that enables multiple
    TimeSeries metadata to be concatenated in an organised fashion.

    Attributes
    ----------
    metadata : `list` of `tuple`
        The list of 3-tuples which each represent a source files metadata.
        The tuples consist of: ( TimeRange, [ colnames ], MetaDict(metadata) )

    Examples
    --------
    >>> from sunpy.timeseries import TimeSeriesMetaData
    >>> from sunpy.time import TimeRange, parse_time
    >>> from sunpy.util import MetaDict
    >>> tr = TimeRange('2012-06-01 00:00','2012-06-02 00:00')
    >>> md = TimeSeriesMetaData(timerange=tr, colnames=['GOES'],
    ...                         meta=MetaDict([('goes_key','goes_val')]))
    >>> tr2 = TimeRange('2012-06-01 12:00','2012-06-02 12:00')
    >>> md.append(tr2, ['EVE'], MetaDict([('eve_key','eve_val')]))
    >>> md.find(parse_time('2012-06-01T21:08:12'))
    |-------------------------------------------------------------------------------------------------|
    |TimeRange                  | Columns         | Meta                                              |
    |-------------------------------------------------------------------------------------------------|
    |2012-06-01 00:00:00        | GOES            | goes_key: goes_val                                |
    |            to             |                 |                                                   |
    |2012-06-02 00:00:00        |                 |                                                   |
    |-------------------------------------------------------------------------------------------------|
    |2012-06-01 12:00:00        | EVE             | eve_key: eve_val                                  |
    |            to             |                 |                                                   |
    |2012-06-02 12:00:00        |                 |                                                   |
    |-------------------------------------------------------------------------------------------------|
    <BLANKLINE>
    >>> md.find(parse_time('2012-06-01T21:08:12')).columns
    ['EVE', 'GOES']
    >>> md.find(parse_time('2012-06-01T21:08:12')).values()
    ['eve_val', 'goes_val']
    >>> md.find(parse_time('2012-06-01T21:08:12')).metas
    [MetaDict([('goes_key', 'goes_val')]), MetaDict([('eve_key', 'eve_val')])]
    >>> md.find(parse_time('2012-06-01T21:08:12'), 'GOES')
    |-------------------------------------------------------------------------------------------------|
    |TimeRange                  | Columns         | Meta                                              |
    |-------------------------------------------------------------------------------------------------|
    |2012-06-01 00:00:00        | GOES            | goes_key: goes_val                                |
    |            to             |                 |                                                   |
    |2012-06-02 00:00:00        |                 |                                                   |
    |-------------------------------------------------------------------------------------------------|

    """

    def __init__(self, meta=None, timerange=None, colnames=None, **kwargs):
        self.metadata = []
        # Parse in arguments
        if not isinstance(meta, type(None)):
            if isinstance(meta, (dict, MetaDict)) and isinstance(timerange, TimeRange) and isinstance(colnames, list):
                # Given a single metadata entry as a dictionary with additional timerange and colnames.
                self.metadata.append((timerange, colnames, meta))
            elif isinstance(meta, tuple):
                # Given a single metadata entry as a tuple.
                self.metadata.append(meta)
            elif isinstance(meta, list):
                # Given a complex metadata list (of tuples)
                self.metadata = copy.copy(meta)
        else:
            # In the event no metadata dictionary is sent we default to something usable
            if isinstance(timerange, TimeRange) and isinstance(colnames, list):
                self.metadata.append((timerange, colnames, MetaDict()))
            elif isinstance(timerange, TimeRange):
                self.metadata.append((timerange, [], MetaDict()))
                warnings.warn("No time range given for metadata. This will mean the metadata can't be linked to columns in data.", Warning)
            else:
                raise ValueError("You cannot create a TimeSeriesMetaData object without specifying a TimeRange")

    def __eq__(self, other):
        """
        Check two TimeSeriesMetaData objects are the same, they have the same
        entries in the same order.

        Parameters
        ----------
        other : `~sunpy.timeseries.metadata.TimeSeriesMetaData`
            The second TimeSeriesMetaData object to compare with.

        Returns
        -------
        result : `bool`
        """
        match = True
        if len(self.metadata) == len(other.metadata):
            for i in range(0,len(self.metadata)):
                #
                if self.metadata[i] != other.metadata[i]:
                    match = False
        else:
            match = False
        return match

    def __ne__(self, other):
        """
        Check two TimeSeriesMetaData objects are not the same, they don't have
        same entries in the same order.

        Parameters
        ----------
        other : `~sunpy.timeseries.GenericTimeSeries`
            The second TimeSeries object to compare with.

        Returns
        -------
        result : `bool`
        """
        return not self == other

    def append(self, timerange, columns, metadata, **kwargs):
        """
        Add the given metadata MetaDict into the metadata list as a tuple with
        it's TimeRange and colnames (list).
        Will add the new entry so the list is in chronological order for the
        TimeRange.start datetime values.

        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            The timerange for which a given metadict is relevant. This will
            generally initilly be the full range of the original file, but if
            the TimeSeries gets truncated this may change appropriately.

        columns : `str`
            A list of the colomn name strings that the metadata is relevant for.

        metadata : `~sunpy.util.metadata.MetaDict` or `OrderedDict` or `dict`
            The dictionary holding the metadata.
        """
        # Parameters
        metadata = MetaDict(metadata)

        # Check the types are correct.
        pos = 0
        if isinstance(timerange, TimeRange):
            for i, meta in enumerate(self.metadata):
                if timerange.start > meta[0].start:
                    pos = i + 1
        else:
            raise ValueError(
                'Incorrect datatime or data for append to TimeSeriesMetaData.')

        # Prepare tuple to append.
        new_metadata = (timerange, columns, metadata)

        # Check this isn't a duplicate entry (same TR and comnames)
        duplicate = False
        if pos < len(self.metadata):
            old_metadata = self.metadata[pos]
            if (new_metadata[0] == old_metadata[0]) and (new_metadata[1] == old_metadata[1]):
                duplicate = True

        # Insert into the given position
        if not duplicate:
            self.metadata.insert(pos, new_metadata)

    def find_indices(self, time=None, colname=None, **kwargs):
        """
        Find the indices for all the metadata entries matching the given filters
        for datetime and/or column name.
        Will return all metadata entry indices if no filters are given.

        Parameters
        ----------
        time : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.

        colname : `str` optional
            A string that can be used to narrow results to specific columns.

        indices : `bool` optional
            If True then return a list of indices, not of MetaDict items.
            Used when other methods use the filters for selecting metadata entries.

        Returns
        -------
        list : `list`
            A list of integers that contain all matching metadata.
        """
        # Parameters
        dt = time
        if not dt:
            dt = False
        elif isinstance(dt, str):
            dt = parse_time(dt)

        # Find all results with suitable timerange.
        results_indices = []
        for i, meta in enumerate(self.metadata):
            if dt in meta[0] or not(dt):
                results_indices.append(i)

        # Filter out only those with the correct column.
        results = []
        for i in results_indices:
            if (colname in self.metadata[i][1]) or (not colname):
                results.append(i)

        return results

    def find(self, time=None, colname=None, **kwargs):
        """
        Find all metadata matching the given filters for datetime and/or column name.
        Will return all metadata entries if no filters are given.

        Parameters
        ----------
        time : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.

        colname : `str` optional
            A string that can be used to narrow results to specific columns.

        Returns
        -------
        metadata : `~sunpy.timeseries.metadata.TimeSeriesMetaData`
            A TimeSeriesMetaData that contain all matching metadata entries.
        """
        # Get the indices
        indices = self.find_indices(time=time, colname=colname, **kwargs)

        # Extract the relevant metadata entries
        metadata = []
        for i in indices:
            metadata.append(copy.copy(self.metadata[i]))

        # Return a TimeSeriesMetaData object
        return TimeSeriesMetaData(meta=metadata)

    def get_index(self, index):
        """
        Return the dictionary entry at the given index.

        Parameters
        ----------
        index : `int`
            The integer index of the metadata entry in the list.

        Returns
        -------
        metadata : `~sunpy.util.metadata.MetaDict`
            An ordered Dictionary containing the metadata at the given index.
        """
        return self.metadata[index][2]

    def get(self, keys, time=None, colname=None, **kwargs):
        """
        Return a TimeSeriesMetaData object of all entries matching the time and
        colname filters with the dictionaries containing only the key value pairs
        with the key matching the given input key.

        Parameters
        ----------
        keys : `str`
            The Key/s to be searched in the dictionary.

        time : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.

        colname : `str` optional
            A string that can be used to narrow results to specific columns.

        itemised : `bool` optional
            Option to allow the return of the time ranges and column names
            (as list) that match each given value.

        Returns
        -------
        metadata : `~sunpy.timeseries.metadata.TimeSeriesMetaData`
            A TimeSeriesMetaData that contain all matching metadata entries but
            with only the requested key/value pairs in the MetaDict objects.
        """
        # Make a list of keys if only one is given
        if isinstance(keys, str):
            keys = [ keys ]

        # Find all metadata entries for the given time/colname filters
        full_metadata = self.find(time=time, colname=colname)
        metadata = []

        # Append to metadata only key:value pairs with requested keys
        for i, entry in enumerate(full_metadata.metadata):
            metadict = MetaDict()
            for curkey, value in entry[2].items():
                for key in keys:
                    if curkey.lower() == key.lower():
                        metadict.update({key:value})
            metadata.append((entry[0], entry[1], metadict))

        # Return a TimeSeriesMetaData object
        return TimeSeriesMetaData(meta=metadata)

    def concatenate(self, tsmetadata2, **kwargs):
        """
        Combine the metadata from a TimeSeriesMetaData object with the current
        TimeSeriesMetaData and return as a new TimeSeriesMetaData object.

        Parameters
        ----------
        tsmetadata2 : `~sunpy.timeseries.TimeSeriesMetaData`
            The second TimeSeriesMetaData object.
        """
        # Create a copy of the metadata
        meta = TimeSeriesMetaData(copy.copy(self.metadata))

        # Append each metadata entry from the second TimeSeriesMetaData object
        # to the original TimeSeriesMetaData object.
        for entry in tsmetadata2.metadata:
            meta.append(entry[0], entry[1], entry[2])

        return meta

    def update(self, dictionary, time=None, colname=None, row=None, overwrite=False, **kwargs):
        """
        Make updates to the MetaDict metadata for all matching metadata entries.

        Parameters
        ----------
        dictionary : `dict` or `OrderedDict` or `~sunpy.util.metadata.MetaDict`
            The second TimeSeriesMetaData object.

        time : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            to filter the metadata entries updated.

        colname : `str` optional
            A string that can be used to narrow results to specific columns.

        overwrite : `bool` optional
            Option to define if the user is able to overwrite already present keys.
            Defaults to False, designed to stop users from being able to
            corrupt/damage the metadict values so easily.
        """
        # Find all matching metadata entries
        indices = self.find_indices(time=time, colname=colname, row=row, indices=True)

        # Now update each matching entries
        for i in indices:
            # Seperate keys for new and current pairs
            old_keys = set(dictionary.keys())
            old_keys.intersection_update(set(self.metadata[i][2].keys()))
            new_keys = set(dictionary.keys())
            new_keys.difference_update(old_keys)

            # Old keys only overwritten if allowed
            for key in (self.metadata[i][2].keys()):
                if key in old_keys and overwrite:
                    self.metadata[i][2][key] = dictionary[key]
            for key in dictionary:
                if key in new_keys:
                    self.metadata[i][2][key] = dictionary[key]

    def _truncate(self, timerange):
        """Removes metadata entries outside of the new (truncated) TimeRange.
        Also adjusts start and end times of time ranges going outside of the
        truncated time range.

        Parameters
        ----------
        timerange : `sunpy.time.TimeRange`
            Either a time range to truncate to.
        """
        truncated = []
        for metatuple in self.metadata:
            # Get metadata time range parameters
            start = metatuple[0].start
            end   = metatuple[0].end
            out_of_range = False

            # Find truncations
            if start < timerange.start and end > timerange.start:
                # Truncate the start
                start = timerange.start
            elif start > timerange.end:
                # Metadata time range starts after truncated data ends.
                out_of_range = True
            if end > timerange.end and start < timerange.end:
                # Truncate the end
                end = timerange.end
            elif end < timerange.start:
                # Metadata time range finishes before truncated data starts.
                out_of_range = True

            # Add the values if applicable
            if not out_of_range:
                truncated.append((TimeRange(start, end), metatuple[1], metatuple[2]))

        # Update the original list
        self.metadata = truncated

    @property
    def columns(self):
        """Returns a list of all the names of the columns in the metadata."""
        all_cols = set()
        for metatuple in self.metadata:
            all_cols.update(metatuple[1])
        all_cols = list(all_cols)
        all_cols.sort()
        return all_cols

    @property
    def metas(self):
        """Returns a list of all the metadict objects in the TimeSeriesMetaData object."""
        all_metas = []
        for metatuple in self.metadata:
            all_metas.append(metatuple[2])
        return all_metas

    @property
    def timeranges(self):
        """Returns a list of all the TimeRange objects the TimeSeriesMetaData object."""
        all_tr = []
        for metatuple in self.metadata:
            all_tr.append(metatuple[0])
        return all_tr

    def values(self):
        """Returns a list of all the values from the metadict objects in each
        entry in the TimeSeriesMetaData object."""
        all_vals = set()
        for metatuple in self.metadata:
            for key, value in metatuple[2].items():
                all_vals.add(str(value))
        all_vals = list(all_vals)
        all_vals.sort()
        return all_vals

    @property
    def time_range(self):
        """Returns the TimeRange of the entire time series meta data."""
        start = self.metadata[0][0].start
        end = self.metadata[0][0].end
        for metatuple in self.metadata:
            if end < metatuple[0].end:
               end = metatuple[0].end
        return TimeRange(start, end)

    def _remove_columns(self, colnames):
        """Removes the given column/s from the TimeSeriesMetaData object.

        Parameters
        ----------
        colnames : `str` or `list`
            The name or names of the columns to be removed.
        """
        # Parameters
        if isinstance(colnames, str):
            colnames = [ colnames ]

        # Create a new list with all metadata entries without colnames
        reduced = []
        for metatuple in self.metadata:
            # Check each colname
            for colname in colnames:
                if colname in metatuple[1]:
                    # Removed from the list.
                    metatuple[1].remove(colname)
            # Add the column if it still has some columns listed
            if len(metatuple[1]) > 0:
                reduced.append(metatuple)

        # Update the original list
        self.metadata = reduced



    def _rename_column(self, old, new):
        """
        Change the name of a column in all the metadata entries.

        Parameters
        ----------
        old : `str`
            The original column name to be changed.

        new : `str`
            The new column name.
        """
        for i in range(0, len(self.metadata)):
            # Update the colnames
            colnames = self.metadata[i][1]
            colnames = [w.replace(old, new) for w in colnames]

            # Replace values
            self.metadata[i] = ( self.metadata[i][0], colnames, self.metadata[i][2] )

    def _validate_meta(self, meta):
        """
        Validate a meta argument.
        """
        # Checking for metadata that may overlap.
        indices = range(0, len(self.metadata))
        for i, j in itertools.combinations(indices, 2):
            # Check if the TimeRanges overlap
            if not ((self.metadata[i][0].end <= self.metadata[j][0].start) or (self.metadata[i][0].start >= self.metadata[j][0].end)):
                # Check column headings overlap
                col_overlap = list(set(self.metadata[i][1]) & set(self.metadata[j][1]))
                # If we have an overlap then show a warning
                if col_overlap:
                    warnings.warn_explicit('Metadata entries ' + str(i) + ' and ' + str(j) + ' contain interleaved data.',
                                           Warning, __file__, inspect.currentframe().f_back.f_lineno)

        # ToDo: Check all entries are in tr.start time order.

        return True

    def to_string(self, depth=10, width=99):
        """
        Print a table-like representation of the TimeSeriesMetaData object.

        Parameters
        ----------
        depth : `int`
            The maximum number of lines to show for each entry. Metadata
            dictionaries and column lists will be truncated if this is small.

        width : `int`
            The number of characters wide to make the entire table.
        """
        # Parameters
        colspace = ' | '
        liswidths = (26, 15, width-2-2*len(colspace) - 26 - 15)
        colheadings = '|' + 'TimeRange'.ljust(100)[:liswidths[0]] + colspace + 'Columns'.ljust(100)[:liswidths[1]] + colspace + 'Meta'.ljust(100)[:liswidths[2]]  + '|'
        rowspace = "-" * (liswidths[0] + len(colspace) + liswidths[1] + len(colspace) + liswidths[2])
        rowspace = '|' + rowspace + '|'

        # Headings
        full = rowspace + '\n' + colheadings + '\n' + rowspace + '\n'

        # Add metadata entries
        for entry in self.metadata:
            # Make lists for each of the columns for each metadata entry
            # Padded to the widths given in liswidths
            lis_range = [ str(entry[0].start), '            to            ', str(entry[0].end) ]
            # Shorten TimeRange representation if depth of only 2
            if depth == 2:
                lis_range = [ str(entry[0].start), str(entry[0].end) ]
            liscols = []
            for col in entry[1]:
                liscols.append(col.ljust(100)[:liswidths[1]])
            lismeta = []
            for key in list(entry[2].keys()):
                string = str(key) + ': ' + str(entry[2][key])
                lismeta.append(string.ljust(100)[:liswidths[2]])

            # Add lines of the entry upto the given depth
            for i in range(0, depth):
                # What to do in the event any of the lists have more entries
                # then the current depth
                if len(lis_range) > i or len(entry[1]) > i or len(lismeta) > i :
                    # The start of the line Str is just a vertical bar/pipe
                    line = '|'
                    # Check we have a time range entry to print
                    if len(lis_range) > i:
                        # Simply add that time range entry to the line Str
                        line += lis_range[i].ljust(100)[:liswidths[0]]
                    else:
                        # No entry to add, so just add a blank space
                        line += ''.ljust(100)[:liswidths[0]]
                    # Add a column break vertical bar/pipe
                    line += colspace
                    # Check we have another column name entry to print
                    if len(entry[1]) > i:
                        # Simply add that column name to the line Str
                        line += entry[1][i].ljust(100)[:liswidths[1]]
                    else:
                        # No entry to add, so just add a blank space
                        line += ''.ljust(100)[:liswidths[1]]
                    # Add a column break vertical bar/pipe
                    line += colspace
                    # Check we have another meta key/value pair to print
                    if len(lismeta) > i:
                        # Simply add that key/value pair to the line Str
                        line += lismeta[i].ljust(100)[:liswidths[2]]
                    else:
                        # No entry to add, so just add a blank space
                        line += ''.ljust(100)[:liswidths[2]]
                    # Finish the line Str with vertical bar/pipe and \n
                    full += line + '|\n'
            # Reached the depth limit, add line to show if the columns are truncated
            if len(lis_range) >= depth or len(entry[1]) >= depth or len(lismeta) >= depth:
                # The start of the line Str is just a vertical bar/pipe
                line = '|'
                # Check we have more time range entries to print
                if len(lis_range) > depth:
                    # We have more time range entries, use ellipsis to show this
                    line += '...'.ljust(100)[:liswidths[0]]
                else:
                    # No entry to add, so just add a blank space
                    line += ''.ljust(100)[:liswidths[0]]
                # Add a column break vertical bar/pipe
                line += colspace
                # Check we have more than one column name entry to print
                if len(entry[1]) > depth:
                    # We have more column name entries, use ellipsis
                    line += '...'.ljust(100)[:liswidths[1]]
                else:
                    # No more column name entries, so just add a blank space
                    line += ''.ljust(100)[:liswidths[1]]
                # Add a column break vertical bar/pipe
                line += colspace
                # Check we have more meta key/value pairs to print
                if len(lismeta) > depth:
                    # We have more key/value pairs, use ellipsis to show this
                    line += '...'.ljust(100)[:liswidths[2]]
                else:
                    # No morekey/value pairs, add a blank space
                    line += ''.ljust(100)[:liswidths[2]]
                # Finish the line Str with vertical bar/pipe and \n
                full += line + '|\n'
            # Add a line to close the table
            full += rowspace + '\n'
        return full

    def __repr__(self):
        return self.to_string()
    def __str__(self):
        return self.to_string()
