# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
__authors__ = ["Alex Hamilton, Stuart Mumford"]
__email__ = "stuart@mumford.me.uk"

from collections import OrderedDict

from sunpy.time import TimeRange, parse_time

import warnings
import inspect

class TimeSeriesMetaData:
    """
    An object used to store metadata for TimeSeries objects that enables multiple
    TimeSeries metadata to be concatenated in an organised fashion.

    Attributes
    ----------
    metadata : `list` of `tuple`
        The list of 3-tuples which each represent a source files metadata.
        The tuples consist of: ( TimeRange, [ colnames ], OrderedDict(metadata) )

    Examples
    --------
    >>> from sunpy.timeseries import TimeSeriesMetaData
    >>> from sunpy.time import TimeRange, parse_time
    >>> from collections import OrderedDict
    >>> tr = TimeRange('2012-06-01 00:00','2012-06-02 00:00')
    >>> md = TimeSeriesMetaData()
    >>> md.append(tr, ['GOES'], OrderedDict([('tr','tr')]))
    >>> md.find(parse_time('2012-06-01T21:08:12'))
    >>> md.find(parse_time('2012-06-01T21:08:12'), 'GOES')   # doctest: +SKIP
    """

    def __init__(self, **kwargs):
        self.metadata = []

    def append(self, timerange, columns, metadata, **kwargs):
        """
        Add the given metadata OrderedDict into the metadata list as a tuple.
        Will add the new entry so the list is in chronological order for the
        TimeRange.start datetime values.
        
        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.
        
        columns : `str`
            A string that can be used to narrow results to specific columns.

        metadata : `OrderedDict`
            A string that can be used to narrow results to specific columns.
        """
        # Check the types are correct.
        pos = len(self.metadata)
        if isinstance(timerange, TimeRange) and isinstance(metadata, OrderedDict):
            for i in range(0, len(self.metadata)):
                if timerange.start < self.metadata[i][0].start:
                    pos = i - 1
        else:
            warnings.warn_explicit("Incorrect datatime or data for append to TimeSeriesMetaData.",
                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
        
        # Prepare tuple to append.
        new_metadata = (timerange, columns, metadata)
        
        # Check this isn't a duplicate entry (same TR and comnames)
        duplicate = False
        if pos < len(self.metadata):
            old_metadata = self.metadata[pos]
            if (new_metadata[0].start == old_metadata[0].start) and (new_metadata[0].end == old_metadata[0].end) and (new_metadata[1] == old_metadata[1]):
                duplicate = True
                
        # Insert into the given position
        if not duplicate:
            self.metadata.insert(pos, new_metadata)
            
    def find(self, datetime=None, colname=None, **kwargs):
        """
        Find all metadata matching the given datetime and optionally column name.
        
        Parameters
        ----------
        datetime : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.
        
        colname : `str` optional
            A string that can be used to narrow results to specific columns.

        indexes : `bool` optional
            If True then return a list of indexes, not of OrderedDict items.
    
        Returns
        -------
        list : `list`
            A list of OrderedDict objects that contain all matching metadata.
        """
        # Parameters
        indexes = kwargs.get('indexes', False)
        dt = datetime
        if not dt:
            dt = False
        elif isinstance(dt, str):
            dt = parse_time(dt)
            
        # Find all results with suitable timerange.
        results_indexes = []
        for i in range(0, len(self.metadata)):
            if dt in self.metadata[i][0] or not(dt):
                results_indexes.append(i)
        
        # Filter out only those with the correct column.
        results = []
        for i in results_indexes:
            if (colname in self.metadata[i][1]) or (not colname):
                if indexes:
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
        metadata : `OrderedDict`
            An ordered Dictionary containing the metadata at the given index.
        """
        return self.metadata[index][1]

    def get(self, key, default=None, datetime=None, colname=None, **kwargs):
        """
        Return a list of all the matching entries in the metadata.
        
        Parameters
        ----------
        key : `str`
            The Key to be searched in the dictionary.

        default : 
            The Value to be returned in case key does not exist.

        datetime : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.
        
        colname : `str` optional
            A string that can be used to narrow results to specific columns.
            
        Returns
        -------
        list : `list`
            Returns a list of the matching entries or the default value.
        """
        # Find all matching metadata entries
        indexes = self.find(datetime, colname, indexes=True)
        
        # Now find each matching entry
        results = []
        for i in indexes:
            # Get the value for the entry in this metadata entry
            value = self.metadata[i][2].get(key)
            # Add to the list if a result was returned
            if value != None:
                results.append(value)
        
        # Return all the results
        return results
        
    def concatenate(self, tsmetadata2, **kwargs):
        """
        Combine the metadata from a TimeSeriesMetaData object with the current
        TimeSeriesMetaData.

        Parameters
        ----------
        tsmetadata2 : `~sunpy.timeseries.TimeSeriesMetaData`
            The second TimeSeriesMetaData object.
        """
        # Append each metadata entry from the second TimeSeriesMetaData object
        # to the original TimeSeriesMetaData object.
        for tuple in tsmetadata2.metadata:
            self.append(tuple[0], tuple[1], tuple[2])

    def update(self, dictionary, datetime=None, colname=None, **kwargs):
        """
        Make updates to the OrderedDict metadata for all matching metadata entries.

        Parameters
        ----------
        dictionary : `dict` or `OrderedDict`
            The second TimeSeriesMetaData object.
            
        datetime : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            to filter the metadata entries updated.
            
        colname : `str` optional
            A string that can be used to narrow results to specific columns.
        """
        # Find all matching metadata entries
        indexes = self.find(datetime, colname, indexes=True)
        
        # Now update each matching entry
        for i in indexes:
            self.metadata[i][2].update(dictionary)

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
        return True
        
if __name__ == "__main__":
    tr_1 = TimeRange('2012-06-01 00:00','2012-06-02 00:00')
    tr_1a = TimeRange('2012-06-01 00:00','2012-06-02 00:00')
    tr_2 = TimeRange('2012-06-02 00:00','2012-06-03 00:00')
    tr_3 = TimeRange('2012-06-03 00:00','2012-06-04 00:00')
    tr_4 = TimeRange('2012-06-04 00:00','2012-06-05 00:00')
    tr_5 = TimeRange('2012-06-02 00:00','2012-06-06 00:00')
    
    # Build a TimeSeriesMetaData object
    md = TimeSeriesMetaData()
    md.append(tr_1, ['GOES'], OrderedDict([('tr_1','tr_1'), ('date-obs', 'yyyy-mm-dd hh-mm')]))
    md.append(tr_2, ['GOES'], OrderedDict([('tr_2','tr_2'), ('date-obs', 'yyyy-mm-dd hh-mm')]))
    md.append(tr_3, ['GOES'], OrderedDict([('tr_3','tr_3'), ('date-obs', 'yyyy-mm-dd hh-mm')]))
    md.append(tr_4, ['GOES'], OrderedDict([('tr_4','tr_4'), ('date-obs', 'yyyy-mm-dd hh-mm')]))
    md.append(tr_5, ['Other'], OrderedDict([('tr_5','tr_5'), ('date-obs', 'yyyy-mm-dd hh-mm')]))

    
    # Check it
    time_1 = parse_time('2012-06-01T21:08:12')
    time_2 = parse_time('2012-06-02T21:08:12')
    time_3 = parse_time('2012-06-03T21:08:12')
    time_4 = parse_time('2012-06-04T21:08:12')
    time_5 = parse_time('2012-06-05T21:08:12')
    time_6 = parse_time('2012-05-05T21:08:12') # Too early
    time_7 = parse_time('2012-06-08T21:08:12') # Too late
    # Without columns
    md.find(time_1) # One result (tr_1)
    md.find(time_2) # Two results (tr_2 and tr_5)
    md.find(time_3) # Two results (tr_5 and tr_3)
    md.find(time_4) # Two results (tr_5 and tr_4)
    md.find(time_5) # One result (tr_5)
    md.find(time_6) # No results (too early)
    md.find(time_7) # No results (too late)
    # With columns
    md.find(time_1, 'Other') # No results (wrong column)
    md.find(time_2, 'Other') # One result (tr_5)
    md.find(time_3, 'Other') # One result (tr_5)
    md.find(time_4, 'Other') # One result (tr_5)
    md.find(time_5, 'Other') # One result (tr_5)
    md.find(time_6, 'Other') # No results (too early)
    md.find(time_7, 'Other')
    md.find(time_1, 'GOES')
    md.find(time_2, 'GOES')
    md.find(time_3, 'GOES')
    md.find(time_4, 'GOES')
    md.find(time_5, 'GOES')
    md.find(time_6, 'GOES') # No results (too early)
    md.find(time_7, 'GOES') # No results (too late)
    
    # renaming columns
    md.rename_column('Other', 'changed')
    md.rename_column('GOES', 'goes')
    
    # Get from the metadata
    date_obs = md.get('date-obs', [ ])
    date_obs = md.get('date-obs', [ ], time_2)
    date_obs = md.get('date-obs', [ ], time_2, 'goes')