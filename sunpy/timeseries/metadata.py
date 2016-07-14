# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
__authors__ = ["Alex Hamilton, Stuart Mumford"]
__email__ = "stuart@mumford.me.uk"

from sunpy.util.metadata import MetaDict
import itertools

import warnings
import inspect

from sunpy.time import TimeRange, parse_time

class TimeSeriesMetaData:
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
    >>> from collections import MetaDict
    >>> tr = TimeRange('2012-06-01 00:00','2012-06-02 00:00')
    >>> md = TimeSeriesMetaData()
    >>> md.append(tr, ['GOES'], MetaDict([('tr','tr')]))
    >>> md.find(parse_time('2012-06-01T21:08:12'))
    >>> md.find(parse_time('2012-06-01T21:08:12'), 'GOES')   # doctest: +SKIP
    """

    def __init__(self, meta=None, timerange=None, colnames=None, **kwargs):
        self.metadata = []
        if meta:
            if isinstance(meta, (dict, MetaDict)) and isinstance(timerange, TimeRange) and isinstance(colnames, list):
                # Given a single metadata entry as a dictionary with additional timerange and colnames.
                self.metadata.append((timerange, colnames, meta))
            elif isinstance(meta, tuple):
                # Given a single metadata entry as a tuple.
                self.metadata.append(meta)
            elif isinstance(meta, list):
                # Given a complex metadata list (of tuples)
                self.metadata = meta.copy()

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
        pos = len(self.metadata)
        if isinstance(timerange, TimeRange):
            for i, meta in enumerate(self.metadata):
                if timerange.start < meta[0].start:
                    pos = i - 1
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
            
    def find(self, time=None, colname=None, row=None, indices=False, **kwargs):
        """
        Find all metadata matching the given filters for datetime and column name.
        Will return all metadata entries if no filters are given.
        
        Parameters
        ----------
        time : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.
        
        row : `int` optional
            Integer index of the row within the data (dataframe) to get the
            datetime index from.
        
        colname : `str` optional
            A string that can be used to narrow results to specific columns.

        indices : `bool` optional
            If True then return a list of indices, not of MetaDict items.
            Used when other methods use the filters for selecting metadata entries.
    
        Returns
        -------
        list : `list`
            A list of MetaDict objects that contain all matching metadata.
        """
        # Parameters
        dt = time
        if row:
            # Given a row index, get the relevant dt.
            dt = self.data.index[row]
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
                if indices:
                    results.append(i)
                else:
                    results.append(self.metadata[i][2])
        
        return results
    
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
        return self.metadata[index][1]

    def get(self, key, default=None, time=None, colname=None, row=None, itemised=False, **kwargs):
        """
        Return a list of all the matching entries in the metadata.
        
        Parameters
        ----------
        key : `str`
            The Key to be searched in the dictionary.

        default : optional
            The Value to be returned in case key does not exist.

        time : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.

        row : `int` optional
            Integer index of the row within the data (dataframe) to get the
            datetime index from.
            
        colname : `str` optional
            A string that can be used to narrow results to specific columns.

        itemised : `bool` optional
            Option to allow the return of the time ranges and column names
            (as list) that match each given value.
            
        Returns
        -------
        list : `list`
            Returns a list of the matching entries or the default value.
        """        
        # Find all matching metadata entries
        indices = self.find(time=time, colname=colname, row=row, indices=True)
        
        # Now find each matching entry
        results = []
        for i in indices:
            # Get the value for the entry in this metadata entry
            value = self.metadata[i][2].get(key)            
            
            # Need more details if itemised.
            if itemised:
                tr = self.metadata[i][0]
                colnames = self.metadata[i][1]
                
                # Append a tuple of these values
                results.append(( tr, colnames, value ))
            else:
                # Add to the list if a result was returned
                if value:
                    results.append(value)
        
        # Remove duplicates.
        if not itemised:
            results = list(set(results))
        
        # Return all the results
        return results
        
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
        meta = TimeSeriesMetaData(self.metadata.copy())
        
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
        indices = self.find(time=time, colname=colname, row=row, indices=True)
        
        # Now update each matching entries
        for i in indices:
            # Should we allow the user to overwrite values?
            if overwrite:
                self.metadata[i][2].update(dictionary)
            else:
                #ToDo: if any new.keys in... iterate through keys, if key is false throw an error otherwise update.
                # Overwrite the original dict over the new to keeps old values.
                old_meta = self.metadata[i][2].copy()
                new_meta = MetaDict(dictionary).copy()
                new_meta.update(old_meta)
                # Now recreate the tuple
                self.metadata[i] = ( self.metadata[i][0], self.metadata[i][1], new_meta )

    def truncate(self, timerange):
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
        return list(all_cols)
        
    def remove_columns(self, colnames):
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



    def rename_column(self, old, new):
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
        #for x, y in itertools.combinations(self.metadata, 2):
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
            listr = [ str(entry[0].start), '            to            ', str(entry[0].end) ]
            liscols = []
            for col in entry[1]:
                liscols.append(col.ljust(100)[:liswidths[1]])
            lismeta = []
            for key in list(entry[2].keys()):
                string = str(key) + ': ' + str(entry[2][key])
                lismeta.append(string.ljust(100)[:liswidths[2]])
            
            # Add lines of the entry upto the given depth
            for i in range(0, depth):
                if len(listr) > i or len(entry[1]) > i or len(lismeta) > i :
                    line = '|'
                    if len(listr) > i:
                        line += listr[i].ljust(100)[:liswidths[0]]
                    else:
                        line += ''.ljust(100)[:liswidths[0]]
                    line += colspace
                    if len(entry[1]) > i:
                        line += entry[1][i].ljust(100)[:liswidths[1]]
                    else:
                        line += ''.ljust(100)[:liswidths[1]]
                    line += colspace
                    if len(lismeta) > i:
                        line += lismeta[i].ljust(100)[:liswidths[2]]
                    else:
                        line += ''.ljust(100)[:liswidths[2]]
                    full += line + '|\n'
            # Add a line to show if the columns are truncated.
            if len(listr) >= depth or len(entry[1]) >= depth or len(lismeta) >= depth:
                line = '|'
                if len(listr) >= depth:
                    line += '...'.ljust(100)[:liswidths[0]]
                else:
                    line += ''.ljust(100)[:liswidths[0]]
                line += colspace
                if len(entry[1]) >= depth:
                    line += '...'.ljust(100)[:liswidths[1]]
                else:
                    line += ''.ljust(100)[:liswidths[1]]
                line += colspace
                if len(lismeta) >= depth:
                    line += '...'.ljust(100)[:liswidths[2]]
                else:
                    line += ''.ljust(100)[:liswidths[2]]
                full += line + '|\n'
            # Add a line to close the table
            full += rowspace + '\n'
        return full

    def __repr__(self):
        return self.to_string()
    def __str__(self):
        return self.to_string()