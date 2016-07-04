# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 18:10:26 2016

@author: alex_
"""

from collections import OrderedDict

from sunpy.time import TimeRange, parse_time

import warnings
import inspect




class TimeSeriesMetaData:
    """
    A metadata object for use with the TimeSeries object.
    """

    def __init__(self, **kwargs):
        self.metadata = []

    def append(self, timerange, columns, metadata, **kwargs):
        """
        Add the given metadata into the metadata table.
        
        Parameters
        ----------
        timerange : `~sunpy.time.TimeRange`
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.
        
        columns : `str`
            A string that can be used to narrow results to specific columns.

        metadata : `OrderedDict`
            A string that can be used to narrow results to specific columns.
    
        Returns
        -------
        list : `list`
            A list of OrderedDict objects that contain all matching metadata.
        """
        # Check the types are correct.
        pos = len(self.metadata)
        if isinstance(timerange, TimeRange) and isinstance(metadata, OrderedDict):
            for i in range(0, len(self.metadata)):
                if timerange.start < self.metadata[i][0].start:
                    pos = i - 1
        else:
            warnings.warn_explicit("Incorrect datatime or data for append to MetaData.",
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
            
    def find(self, datetime, colname=None):
        """
        Find all metadata matching the given date/time and optional column name.
        
        Parameters
        ----------
        datetime : `str` or `~datetime.datetime`
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.
        
        colname : `str` optional
            A string that can be used to narrow results to specific columns.
    
        Returns
        -------
        list : `list`
            A list of OrderedDict objects that contain all matching metadata.
        """
        # Extract the datetime object.
        dt = datetime
        if isinstance(dt, str):
            dt = parse_time(dt)
            
        # Find all results with suitable timerange.
        results_indexes = []
        for i in range(0, len(self.metadata)):
            if dt in self.metadata[i][0]:
                results_indexes.append(i)
        
        # Filter out only those with the correct column.
        results = []
        for i in results_indexes:
            #print('\n\nself.metadata[i][1]: ' + str(self.metadata[i][1]) + '\n\n')
            if (colname in self.metadata[i][1]) or (not colname):
                results.append(self.metadata[i][2])
        
        return results
    
    def get(self, index):
        # Return the dictionary entry at the given index.
        return self.metadata[index][1]
        
    def concatenate(self, tsmetadata2, **kwargs):
        """
        Combine the metadata from a TimeSeriesMetaData object with the current
        metadata.
        """
        for tuple in tsmetadata2.metadata:
            self.append(tuple[0], tuple[1], tuple[2])
        
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
    
    # Build a MetaData object
    md = TimeSeriesMetaData()
    md.append(tr_1, ['GOES'], OrderedDict([('tr_1','tr_1')]))
    md.append(tr_2, ['GOES'], OrderedDict([('tr_2','tr_2')]))
    md.append(tr_3, ['GOES'], OrderedDict([('tr_3','tr_3')]))
    md.append(tr_4, ['GOES'], OrderedDict([('tr_4','tr_4')]))
    md.append(tr_5, ['Other'], OrderedDict([('tr_5','tr_5')]))
    
    #metadict = OrderedDict([(tr_1, 'tr_1'), (tr_2, 'tr_2'), (tr_3, 'tr_3'), (tr_4, 'tr_4'), (tr_5, 'tr_5')])
    #md = MetaData(metadict)
    
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