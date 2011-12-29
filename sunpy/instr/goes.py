# -*- coding: utf-8 -*-
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
#
# <License info will go here...>

from __future__ import absolute_import
import matplotlib.pyplot as plt
import matplotlib.dates
import datetime

"""Provides programs to process and analyze GOES data.

"""
#TODO goes should really be an object which inherits from a timeseries object

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
    ax.set_xlabel(datetime.datetime.isoformat(t[0])[0:10])
    
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