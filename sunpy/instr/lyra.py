from __future__ import absolute_import

import datetime
import urlparse
import calendar
import sqlite3
import sunpy.net.download

def download_lytaf_database(lytaf_dir=''):
    """download the latest version of the Proba-2 pointing database from the Proba2 Science Center"""
    dl=sunpy.net.download.Downloader()
    dl.download('http://proba2.oma.be/lyra/data/lytaf/annotation_ppt.db',path=lytaf_dir)
    
    print 'LYTAF update completed'
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
    conn=sqlite3.connect(lytaf_dir + 'annotation_ppt.db')
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
        lar_entry={'roi_discription':'LYRA LYTAF event',
        'start_time':datetime.datetime.utcfromtimestamp(l[1]),
        'ref_time':datetime.datetime.utcfromtimestamp(l[2]),
        'end_time':datetime.datetime.utcfromtimestamp(l[3]),
        'event_type_id':l[4],'event_type_description':_lytaf_event2string(l[4])}

        #create output tuple for each entry in list
        #entry=(insertion_time,start_time,ref_time,end_time,event_type,event_type_info[0])
        #output a list of dictionaries
        output.append(lar_entry)
    
    return output

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
