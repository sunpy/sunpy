"""
    Provides programs to process and analyze GOES data. This module is
    currently in development.
   
    Examples
    --------
    To make a GOES plot
    >>> import sunpy.instr.goes as goes
    >>> f = goes.get_goes_file(['2011/04/04', '2011/04/05'])
    >>> data = sunpy.instr.goes.parse_goes(f[0])
    >>> goes.show(res.get('time_tag'), res.get('A_FLUX'), res.get('B_FLUX'))
    
    To Do
    -----
    * An object should be developed to encapsulate all of goes functionality
    * It currently looks like the GOES API only allows for downloading an entire
        days worth of data at full cadence. This is a lot of data to work with.
        Plotting it takes too long. Should probably sparsify it.
    
    See Also
    --------
    For current data plots see http://www.swpc.noaa.gov/rt_plots/xray_5m.html
    
    References
    ----------
    | http://www.ngdc.noaa.gov/goes/sem/getData

    """

from __future__ import absolute_import
import matplotlib.pyplot as plt
import matplotlib.dates
from datetime import datetime
from sunpy.time import TimeRange
from sunpy.util.util import break_time
from sunpy.util.util import anytim
import urllib
import csv

def get_goes_file(time_range):
    """Get the goes data"""
    
    _time_range = TimeRange(time_range)
    
    url_root = 'http://www.ngdc.noaa.gov/goes/sem/getData/goes15/xrs_2s.csv?fromDate='
    
    url = url_root + break_time(_time_range.t1)[0:8] + '&toDate=' + break_time(_time_range.t2)[0:8]
    url = url + '&file=true'
    print('Downloading file: ' + url)
    f = urllib.urlretrieve(url)

    return f

def parse_goes(file):
    """Parse a goes file"""
    
    reader = csv.reader(open(file, "rb"), delimiter = ',', skipinitialspace = True)
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
    """Plot the standard GOES plot."""
        
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