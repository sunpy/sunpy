# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
__authors__ = ["Alex Hamilton, Stuart Mumford"]
__email__ = "stuart@mumford.me.uk"

from sunpy.util.metadata import MetaDict

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
            if (new_metadata[0].start == old_metadata[0].start) and (new_metadata[0].end == old_metadata[0].end) and (new_metadata[1] == old_metadata[1]):
                duplicate = True
                
        # Insert into the given position
        if not duplicate:
            self.metadata.insert(pos, new_metadata)
            
    def find(self, datetime=None, colname=None, indexes=False, **kwargs):
        """
        Find all metadata matching the given filters for datetime and column name.
        Will return all metadata entries if no filters are given.
        
        Parameters
        ----------
        datetime : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.
        
        colname : `str` optional
            A string that can be used to narrow results to specific columns.

        indexes : `bool` optional
            If True then return a list of indexes, not of MetaDict items.
            Used when other methods use the filters for selecting metadata entries.
    
        Returns
        -------
        list : `list`
            A list of MetaDict objects that contain all matching metadata.
        """
        # Parameters
        dt = datetime
        if not dt:
            dt = False
        elif isinstance(dt, str):
            dt = parse_time(dt)
            
        # Find all results with suitable timerange.
        results_indexes = []
        for i, meta in enumerate(self.metadata):
            if dt in meta[0] or not(dt):
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
        metadata : `~sunpy.util.metadata.MetaDict`
            An ordered Dictionary containing the metadata at the given index.
        """
        return self.metadata[index][1]

    def get(self, key, default=None, datetime=None, colname=None, itemised=False, **kwargs):
        """
        Return a list of all the matching entries in the metadata.
        
        Parameters
        ----------
        key : `str`
            The Key to be searched in the dictionary.

        default : optional
            The Value to be returned in case key does not exist.

        datetime : `str` or `~datetime.datetime` optional
            The string (parsed using the `~sunpy.time.parse_time`) or datetime
            that you need metadata for.
        
        colname : `str` optional
            A string that can be used to narrow results to specific columns.

        itemised : `bool` optional
            Option to allow the return of the time ranges and column names (as list) that match each given value.
            
        Returns
        -------
        list : `list`
            Returns a list of the matching entries or the default value.
        """        
        # Find all matching metadata entries
        indexes = self.find(datetime=datetime, colname=colname, indexes=True)
        
        # Now find each matching entry
        results = []
        for i in indexes:
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

    def update(self, dictionary, datetime=None, colname=None, overwrite=False, **kwargs):
        """
        Make updates to the MetaDict metadata for all matching metadata entries.

        Parameters
        ----------
        dictionary : `dict` or `OrderedDict` or `~sunpy.util.metadata.MetaDict`
            The second TimeSeriesMetaData object.
            
        datetime : `str` or `~datetime.datetime` optional
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
        indexes = self.find(datetime=datetime, colname=colname, indexes=True)
        
        # Now update each matching entries
        for i in indexes:
            # Should we allow the user to overwrite values?
            if overwrite:
                self.metadata[i][2].update(dictionary)
            else:
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
    md.append(tr_1, ['GOES'], MetaDict([('tr_1','tr_1'), ('date-obs', 'yyyy-mm-dd hh-mm')]))
    md.append(tr_2, ['GOES'], MetaDict([('tr_2','tr_2'), ('date-obs', 'yyyy-mm-dd hh-mm')]))
    md.append(tr_3, ['GOES'], MetaDict([('tr_3','tr_3'), ('date-obs', 'yyyy-mm-dd hh-mm')]))
    md.append(tr_4, ['GOES'], MetaDict([('tr_4','tr_4'), ('date-obs', 'yyyy-mm-dd hh-mm')]))
    md.append(tr_5, ['Other'], MetaDict([('tr_5','tr_5'), ('date-obs', 'yyyy-mm-dd hh-mm')]))
    
    # Check it
    time_1 = parse_time('2012-06-01T21:08:12')
    time_2 = parse_time('2012-06-02T21:08:12')
    time_3 = parse_time('2012-06-03T21:08:12')
    time_4 = parse_time('2012-06-04T21:08:12')
    time_5 = parse_time('2012-06-05T21:08:12')
    time_6 = parse_time('2012-05-05T21:08:12') # Too early
    time_7 = parse_time('2012-06-08T21:08:12') # Too late
    # Without columns
    md.find(datetime=time_1) # One result (tr_1)
    md.find(datetime=time_2) # Two results (tr_2 and tr_5)
    md.find(datetime=time_3) # Two results (tr_5 and tr_3)
    md.find(datetime=time_4) # Two results (tr_5 and tr_4)
    md.find(datetime=time_5) # One result (tr_5)
    md.find(datetime=time_6) # No results (too early)
    md.find(datetime=time_7) # No results (too late)
    # With columns
    md.find(datetime=time_1, colname='Other') # No results (wrong column)
    md.find(datetime=time_2, colname='Other') # One result (tr_5)
    md.find(datetime=time_3, colname='Other') # One result (tr_5)
    md.find(datetime=time_4, colname='Other') # One result (tr_5)
    md.find(datetime=time_5, colname='Other') # One result (tr_5)
    md.find(datetime=time_6, colname='Other') # No results (too early)
    md.find(datetime=time_7, colname='Other')
    md.find(datetime=time_1, colname='GOES')
    md.find(datetime=time_2, colname='GOES')
    md.find(datetime=time_3, colname='GOES')
    md.find(datetime=time_4, colname='GOES')
    md.find(datetime=time_5, colname='GOES')
    md.find(datetime=time_6, colname='GOES') # No results (too early)
    md.find(datetime=time_7, colname='GOES') # No results (too late)
    
    # Get from the metadata
    date_obs = md.get('date-obs', [ ])
    date_obs = md.get('date-obs', [ ], datetime=time_2)
    date_obs = md.get('date-obs', [ ], datetime=time_2, colname='GOES')
    date_obs = md.get('date-obs', itemised=True)
    
    # Update
    md.update({'new_key_1':'added to all.'})
    md.update({'new_key_2':'added to all at time_3.'}, datetime=time_3)
    md.update({'new_key_3':'added only to time_4'}, datetime=time_4, colname=None)
    md.update({'comment':'short'}, overwrite=True)

    # renaming columns
    md.rename_column('Other', 'changed')
    md.rename_column('changed', 'Other')
    md.rename_column('GOES', 'goes')
    md.rename_column('goes', 'GOES')