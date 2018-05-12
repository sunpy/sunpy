"""
================================================
Parsing times with sunpy.time.parse_time
================================================

Example to show some example usage of parse_time
"""

from datetime import datetime, date
import numpy as np
import pandas

from sunpy.time import parse_time

# dict used for coloring the terminal output
col = {'y': '\x1b[93m', 'g': '\x1b[92m', 'r': '\x1b[96m', 'bold': '\x1b[1m',
       'end': '\x1b[0m'}


def print_time(*args, **kwargs):
    '''Parses and pretty prints a parse_time compatible object
    '''

    # Parse the time
    time = parse_time(*args, **kwargs)  # Pass all arguments to parse_time

    # Color and print to terminal
    print(col['r'] + '\nInput string/object: ' + col['end'] +
          col['bold'] + '{ts!r}'.format(ts=args[0])+col['end'])
    print(col['r'] + 'Parsed Time: ' + col['end'] + col['y'] + col['bold'] +
          '{time!r}'.format(time=time) + col['end'])


# Strings
print('\nSTRINGS')
print_time('2005-08-04T00:18:02.000', scale='tt')
print_time('20140101000001')
print_time('2016.05.04_21:08:12_TAI')
print_time('1995-12-31 23:59:60')  # Leap second
print_time('1995-Dec-31 23:59:60')

# datetime
print('\nDATETIME')
print_time(datetime.now(), scale='tai')
print_time(date.today())

# numpy
print('\nnumpy.datetime64')
print_time(np.datetime64('1995-12-31 18:59:59-0500'))
print_time(np.arange('2005-02-01T00', '2005-02-01T10', dtype='datetime64'))

# astropy compatible times
print('\nAstroPy compatible')
print_time(1234.0, format='jd')
print_time('B1950.0', format='byear_str')
print_time('2001-03-22 00:01:44.732327132980', scale='utc',
           location=('120d', '40d'))  # pass location

# pandas
print_time(pandas.Timestamp(datetime(1966, 2, 3)))
print_time(
    pandas.Series([[datetime(2012, 1, 1, 0, 0),
                    datetime(2012, 1, 2, 0, 0)],
                   [datetime(2012, 1, 3, 0, 0),
                    datetime(2012, 1, 4, 0, 0)]]))
