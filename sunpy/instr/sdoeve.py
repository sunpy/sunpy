# -*- coding: utf-8 -*-
#
# Author: Steven Christe <steven.d.christe@nasa.gov>
#
# <License info will go here...>

from __future__ import absolute_import
import urllib
import csv
from datetime import datetime, date, time  
import sunpy.util

"""
    Provides programs to process and analyze EVE data.

"""  

def get_latest_l0cs_goes_data():
    """Grab the latest EVE GOES Proxy data and plot it in a standard 
    (GOES) plot format"""
    #TODO should this be in the net module?
    
    url = 'http://lasp.colorado.edu/eve/data_access/quicklook/quicklook_data/L0CS/LATEST_EVE_L0CS_DIODES_1m.txt'
    
    f = urllib.urlretrieve(url)
    reader = csv.reader(open(f[0], "rb"), delimiter = ' ', skipinitialspace = True)
    
    i = 0
    
    t = []
    xrsb = []
    xrsa = []
    
    for row in reader:
        if row[0][0] != ';':
            #read the date line
            if i == 0:
                d = date(int(row[0]),int(row[2]),int(row[3]))
            else:
                 t.append(time(int(row[0][0:2]),int(row[0][2:4])))
                 xrsb.append(float(row[1]))
                 xrsa.append(float(row[2])) 
            i = i + 1
           
    ts = [datetime.combine(d,s) for s in t]
    
    return [ts,xrsa, xrsb]

def show_latest_l0cs_goes_data():
    """Plot the latest EVE GOES proxy data in a standard GOES plot."""
    from sunpy.instr.goes import show as goes_show
    
    data = get_latest_l0cs_goes_data()
    
    goes_show(data[0], data[1], data[2], 
              title = 'EVE GOES Proxy Xray Flux (1 minute data)')
    
def get_l0cs_data(time_range):
    return 0
    
def get_l0cs_date(date, kind = None):
    from sunpy.util.util import anytim
    
    url_root = 'http://lasp.colorado.edu/eve/data/quicklook/L0CS/SpWx/'
    _date = anytim(date)
    
    url = url_root + _date.strftime('%Y/%Y%m%d') + '_EVE_L0CS_DIODES_1m.txt'
    url_counts = url_root + _date.strftime('%Y/%Y%m%d') + 
        '_EVE_L0CS_DIODES_1m_counts.txt'
    
    f = urllib.urlretrieve(url)
    reader = csv.reader(open(f[0], "rb"), delimiter = ' ', skipinitialspace = True)
  
    field_names = ('hhmm', 'xrs-b', 'xrs-a', 'sem', 'ESPquad', 'esp171', 
                   'esp257', 'esp304', 'esp366', 'espdark', 'megsp', 'megsdark', 
                   'q0esp', 'q1esp', 'q3esp', 'cmlat', 'cmlon')
    
    t = []

    for row in reader:
        if row[0][0] != ';':
            #read the date line
            if i == 0:
                d = date(int(row[0]),int(row[2]),int(row[3]))
            else:
                 t.append(time(int(row[0][0:2]),int(row[0][2:4])))
                 xrsb.append(float(row[1]))
                 xrsa.append(float(row[2])) 
            i = i + 1
   