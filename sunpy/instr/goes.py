from __future__ import absolute_import

import os
from sunpy.net import hek
from sunpy.time import parse_time

__all__ = ['get_goes_event_list']

def get_goes_event_list(trange,goes_class_filter=None):
    """A function to retrieve a list of flares detected by GOES within a given time range.

    Parameters
    ----------
    trange: a SunPy TimeRange object

    goes_class_filter: (optional) string
        a string specifying a minimum GOES class for inclusion in the list, e.g. 'M1', 'X2'."""
    
    
    #use HEK module to search for GOES events
    client=hek.HEKClient()
    event_type='FL'
    tstart=trange.start()
    tend=trange.end()

    #query the HEK for a list of events detected by the GOES instrument between tstart and tend (using a GOES-class filter)
    if goes_class_filter:
        result=client.query(hek.attrs.Time(tstart,tend),hek.attrs.EventType(event_type),hek.attrs.FL.GOESCls > goes_class_filter,hek.attrs.OBS.Observatory == 'GOES')
    else:
        result=client.query(hek.attrs.Time(tstart,tend),hek.attrs.EventType(event_type),hek.attrs.OBS.Observatory == 'GOES')

    #want to condense the results of the query into a more manageable dictionary
    #keep event data, start time, peak time, end time, GOES-class, location, active region source (as per GOES list standard)
    #make this into a list of dictionaries
    goes_event_list=[]

    for r in result:
        goes_event={}
        goes_event['event_date'] = parse_time(r['event_starttime']).date().strftime('%Y-%m-%d')
        goes_event['start_time'] =parse_time(r['event_starttime'])
        goes_event['peak_time'] = parse_time(r['event_peaktime'])
        goes_event['end_time'] = parse_time(r['event_endtime'])
        goes_event['goes_class'] = str(r['fl_goescls'])
        goes_event['goes_location'] = r['event_coord1'],r['event_coord2']
        goes_event['noaa_active_region'] = r['ar_noaanum']
        goes_event_list.append(goes_event)

    return goes_event_list
