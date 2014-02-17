from __future__ import absolute_import

import datetime
import calendar
import sqlite3
import numpy as np
import os
from sunpy.time import parse_time
import urllib

def download_lytaf_database(lytaf_dir=''):
    """download the latest version of the Proba-2 pointing database from the Proba2 Science Center"""
    #dl=sunpy.net.download.Downloader()
    #dl.download('http://proba2.oma.be/lyra/data/lytaf/annotation_ppt.db',path=lytaf_dir)
    url='http://proba2.oma.be/lyra/data/lytaf/annotation_ppt.db'
    destination=os.path.join(lytaf_dir,'annotation_ppt.db')
    urllib.urlretrieve(url,destination)
    
    return


def get_lytaf_events(timerange,lytaf_dir=''):
    """returns a list of LYRA pointing events that occured during a timerange"""
    #timerange is a TimeRange object
    #start_ts and end_ts need to be unix timestamps
    st_timerange = timerange.start()
    start_ts=calendar.timegm(st_timerange.timetuple())
    en_timerange=timerange.end()
    end_ts=calendar.timegm(en_timerange.timetuple())
    
    #involves executing SQLite commands from within python.
    #connect to the SQlite database
    #conn=sqlite3.connect(lytaf_dir + 'annotation_ppt.db')
    conn=sqlite3.connect(os.path.join(lytaf_dir,'annotation_ppt.db'))
    cursor=conn.cursor()

    #create a substitute tuple out of the start and end times for using in the database query
    query_tup=(start_ts,end_ts,start_ts,end_ts,start_ts,end_ts)

    #search only for events within the time range of interest (the lightcurve start/end). Return records ordered by start time
    result=(cursor.execute('select * from event WHERE((begin_time > ? AND begin_time < ?) OR (end_time > ? AND end_time < ?)' +  
                'OR (begin_time < ? AND end_time > ?)) ORDER BY begin_time ASC',query_tup))
    
    #get all records from the query in python list format. 
    list=result.fetchall()

    #times are in unix time - want to use datetime instead
    output=[]

    for l in list:
        #create a dictionary for each entry of interest
        lar_entry={'roi_description':'LYRA LYTAF event',
        'start_time':datetime.datetime.utcfromtimestamp(l[1]),
        'ref_time':datetime.datetime.utcfromtimestamp(l[2]),
        'end_time':datetime.datetime.utcfromtimestamp(l[3]),
        'event_type_id':l[4],'event_type_description':_lytaf_event2string(l[4])[0]}

        #create output tuple for each entry in list
        #entry=(insertion_time,start_time,ref_time,end_time,event_type,event_type_info[0])
        #output a list of dictionaries
        output.append(lar_entry)
    
    return output

def split_series_using_lytaf(timearray,data,lar):
    """
    Proba-2 analysis code for splitting up LYRA timeseries around locations where LARs
    (and other data events) are observed.

    Inputs
    ------
    timearray - an array of times that can be understood by the SunPy parse_time function
    data - data array corresponding to the given time array
    lar - list of events obtained from querying the LYTAF database using lyra.get_lytaf_events()

    Output
    ------
    A list of dictionaries. Each dictionary contains a sub-series corresponding to an interval of 'good data'
    """    
    #lar is a dictionary with tags:
    #'start_time'
    #'end_time'
    #'ref_time'
    #'roi_description'
    #'event_type_description'
    #'event_type_id'

    
    n=len(timearray)
    mask=np.ones(n)
    el=len(lar)

    #make the input time array a list of datetime objects
    datetime_array=[]
    for tim in timearray:
        datetime_array.append(parse_time(tim))
        

        #scan through each entry retrieved from the LYTAF database
    for j in range(0,el):
        #want to mark all times with events as bad in the mask, i.e. = 0
        start_dt=lar[j]['start_time']
        end_dt=lar[j]['end_time']

        #find the start and end indices for each event
        start_ind=np.searchsorted(datetime_array,start_dt)
        end_ind=np.searchsorted(datetime_array,end_dt)

        #append the mask to mark event as 'bad'
        mask[start_ind:end_ind] = 0


    diffmask=np.diff(mask)
    tmp_discontinuity=np.where(diffmask != 0.)
    #disc contains the indices of mask where there are discontinuities
    disc = tmp_discontinuity[0]

    if len(disc) == 0:
        print 'No events found within time series interval. Returning original series.'
        return [{'subtimes':datetime_array,'subdata':data}]
    
    #-1 in diffmask means went from good data to bad
    #+1 means went from bad data to good

    #want to get the data between a +1 and the next -1

    #if the first discontinuity is a -1 then the start of the series was good. 
    if diffmask[disc[0]] == -1.0:
        #make sure we can always start from disc[0] below
        disc=np.insert(disc,0,0)
    
    split_series=[]

    limit=len(disc)
    #now extract the good data regions and ignore the bad ones
    for h in range(0,limit,2):

        if h == limit-1:
            #can't index h+1 here. Go to end of series
            subtimes=datetime_array[disc[h]:-1]
            subdata=data[disc[h]:-1]
            subseries={'subtimes':subtimes,'subdata':subdata}
            split_series.append(subseries)
        else:
            subtimes=datetime_array[disc[h]:disc[h+1]]
            subdata=data[disc[h]:disc[h+1]]
            subseries={'subtimes':subtimes,'subdata':subdata}
            split_series.append(subseries)

    return split_series

def _lytaf_event2string(integers):
    if type(integers) == int:
        integers=[integers]
    #else:
    #    n=len(integers)
    out=[]

    for i in integers:
        if i == 1:
            out.append('LAR')
        if i == 2:
            out.append('N/A')
        if i == 3:
            out.append('UV occult.')
        if i == 4:
            out.append('Vis. occult.')
        if i == 5:
            out.append('Offpoint')
        if i == 6:
            out.append('SAA')
        if i == 7:
            out.append('Auroral zone')
        if i == 8:
            out.append('Moon in LYRA')
        if i == 9:
            out.append('Moon in SWAP')
        if i == 10:
            out.append('Venus in LYRA')
        if i == 11:
            out.append('Venus in SWAP')

    return out
