# -*- coding: utf-8 -*-
#
#
# <License info will go here...>
"""
    Provides programs to process and analyze GOES data. 
    
    .. warning:: This module is still in development!
    
    .. todo:: An object should be developed to encapsulate all of goes functionality
    
    .. todo:: The current GOES API can only provide a full day worth of data.
    
"""

from __future__ import absolute_import
import matplotlib.pyplot as plt
import matplotlib.dates
from datetime import datetime
from sunpy.time import TimeRange
import urllib
import csv

def get_file(time_range):
    """Download the GOES data through the GOES SEM - DataService Web API. 
    
    Parameters
    ----------
    time_range : A TimeRange or time range compatible string

    Returns
    -------
    value : tuple
        Return a tuple (filename, headers) where filename is the local file 
        name under which the object can be found, and headers is 
        whatever the info() method of the object returned by urlopen.

    See Also
    --------

    Examples
    --------
    >>> import sunpy.instr.goes as goes
    >>> goes.get_file(('2011/04/04', '2011/04/05'))
    
    Reference
    ---------
    | http://www.ngdc.noaa.gov/goes/sem/getData
    
    .. note:: This API is currently limited to providing data from 
    whole days only.

    """
    
    _time_range = TimeRange(time_range)
    
    url_root = 'http://www.ngdc.noaa.gov/goes/sem/getData/goes15/xrs_2s.csv?fromDate='
    
    url = url_root + _time_range.t1.strftime("%Y%m%d") + '&toDate=' + _time_range.t2.strftime("%Y%m%d")
    url = url + '&file=true'
    print('Downloading file: ' + url)
    f = urllib.urlretrieve(url)

    return f

def parse_file(filename):
    """Parse a GOES file.
    
    Parameters
    ----------
    filename : The filename of a GOES csv file.

    Returns
    -------
    value : dict
        Returns a dictionary.

    See Also
    --------

    Examples
    --------
    >>> import sunpy.instr.goes as goes
    >>> f = goes.get_file(('2011/04/04', '2011/04/05'))
    >>> goes.parse_file(f[0])

    """
    
    reader = csv.reader(open(filename, "rb"), delimiter = ',', skipinitialspace = True)
    headerline = reader.next()
    
    t = []
    time_tag_ms = []
    a_qual_flag = []
    a_count = []
    a_flux = []
    b_qual_flag = []
    b_count = []
    b_flux = []
    
    for row in reader:
        t.append(datetime.strptime(row[0], '%Y-%m-%d %H:%M:%S.%f'))
        time_tag_ms.append(row[1])
        a_qual_flag.append(row[2])
        a_count.append(row[3])
        a_flux.append(row[4])
        b_qual_flag.append(row[5])
        b_count.append(row[6])
        b_flux.append(row[7])
    
    result = {headerline[0]: t, headerline[1]: time_tag_ms, 
              headerline[2]: a_qual_flag, headerline[3]: a_count, 
              headerline[4]: a_flux, headerline[5]: b_qual_flag, 
              headerline[6]: b_count, headerline[7]: b_flux}
    
    return result

def show(t, xrsa, xrsb, title = 'GOES Xray Flux'):
    """Plot a standard GOES plot.
    
    Parameters
    ----------
    t : tuple
        A tuple containing the times of the measurements.
    xrsa : tuple
        A tuple containing the GOES low energy channel .
    xrsb : tuple
        A tuple containing the GOES high energy channel.
    title : string
        The title of the plot

    Returns
    -------
    value : dict
        Returns a dictionary.

    See Also
    --------

    Examples
    --------
    >>> import sunpy.instr.goes as goes
    >>> f = goes.get_file(('2011/04/04', '2011/04/05'))
    >>> data = goes.parse_file(f[0])
    >>> goes.show(data.get('time_tag'), data.get('A_FLUX'), data.get('B_FLUX'))

    """
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    dates = matplotlib.dates.date2num(t)

    ax.plot_date(dates, xrsa, '-', label = '0.5--4.0 $\AA$', color = 'blue', lw = 2)
    ax.plot_date(dates, xrsb, '-', label = '1.0--8.0 $\AA$', color = 'red', lw = 2)
    ax.set_yscale("log")
    ax.set_ylim(1e-9, 1e-2)
    ax.set_title(title)
    ax.set_ylabel('Watts m$^{-2}$')
    ax.set_xlabel(datetime.isoformat(t[0])[0:10])
    
    ax2 = ax.twinx()
    ax2.set_yscale("log")
    ax2.set_ylim(1e-9, 1e-2)
    ax2.set_yticks((1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
    ax2.set_yticklabels((' ','A','B','C','M','X',' '))
    
    ax.yaxis.grid(True, 'major')
    ax.xaxis.grid(False, 'major')
    ax.legend()
    
    formatter = matplotlib.dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(formatter)
    
    ax.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M')
    fig.autofmt_xdate()
    fig.show()