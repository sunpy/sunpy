# -*- coding: utf-8 -*-
"""Provides programs to process and analyze PROBA2/LYRA data."""
from __future__ import absolute_import

import os
import datetime
import urlparse
import calendar
import sqlite3
import subprocess

from matplotlib import pyplot as plt
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import pandas

import sunpy
from sunpy.lightcurve import LightCurve 
from sunpy.time import parse_time

__all__ = ['LYRALightCurve']

class LYRALightCurve(LightCurve):
    """LYRA light curve definition

    Examples
    --------
    import sunpy
    lyra = sunpy.lightcurve.LYRALightCurve.create()
    lyra = sunpy.lightcurve.LYRALightCurve.create('~/Data/lyra/lyra_20110810-000000_lev2_std.fits')
    lyra = sunpy.lightcurve.LYRALightCurve.create('2011/08/10')
    lyra = sunpy.lightcurve.LYRALightCurve.create("http://proba2.oma.be/lyra/data/bsd/2011/08/10/lyra_20110810-000000_lev2_std.fits")
    lyra.peek()

    References
    ----------
    | http://proba2.sidc.be/data/LYRA
    """
    #lytaf is blank by default. It can be populated with the get_lytaf_as_list method.
    lytaf=None

    def peek(self, names=3, **kwargs):
        """Plots the LYRA data

        See: http://pandas.sourceforge.net/visualization.html
        """
        lyranames = (('Lyman alpha','Herzberg cont.','Al filter','Zr filter'),
                 ('120-123nm','190-222nm','17-80nm + <5nm','6-20nm + <2nm'))

        # Choose title if none was specified
        #if not kwargs.has_key("title"):
        #    if len(self.data.columns) > 1:
        #        kwargs['title'] = 'LYRA data'
        #    else:
        #        if self._filename is not None:
        #            base = self._filename
        #            kwargs['title'] = os.path.splitext(base)[0]
        #        else:
        #            kwargs['title'] = 'LYRA data'

        """Shows a plot of all four light curves"""
        figure = plt.figure()
        axes = plt.gca()

        axes = self.data.plot(ax=axes, subplots=True, sharex=True, **kwargs)
        #plt.legend(loc='best')

        for i, name in enumerate(self.data.columns):
            if names < 3:
                name = lyranames[names][i]
            else:
                name = lyranames[0][i] + ' (' + lyranames[1][i] + ')'
            axes[i].set_ylabel("%s (%s)" % (name, "W/m**2"))

        axes[0].set_title("LYRA ("+ self.data.index[0].strftime('%Y-%m-%d') +")")
        axes[-1].set_xlabel("Time")

        figure.show()

        return figure

    def download_lytaf_database(self,dir=''):
        """download the latest version of the Proba-2 pointing database from the Proba2 Science Center"""
        subprocess.call(['curl','http://proba2.oma.be/lyra/data/lytaf/annotation_ppt.db','-o',dir+'annotation_ppt.db'])

        print 'LYTAF update completed'
        return


    def get_lytaf_events(self,database_dir=''):
        """populates self.lytaf with a list of LYRA pointing events that occured during the lightcurve timerange"""
        #timerange is a TimeRange object
        timerange=self.time_range()
        #start_ts and end_ts need to be unix timestamps
        s=timerange.start()
        start_ts=calendar.timegm(s.timetuple())
        e=timerange.end()
        end_ts=calendar.timegm(e.timetuple())
    
        #involves executing SQLite commands from within python.
        #connect to the SQlite database
        conn=sqlite3.connect(database_dir+'annotation_ppt.db')
        c=conn.cursor()

        #create a substitute tuple out of the start and end times for using in the database query
        tup=(start_ts,end_ts,start_ts,end_ts,start_ts,end_ts)

        #search only for events within the time range of interest (the lightcurve start/end). Return records ordered by start time
        r=(c.execute('select * from event WHERE((begin_time > ? AND begin_time < ?) OR (end_time > ? AND end_time < ?)' +  
                'OR (begin_time < ? AND end_time > ?)) ORDER BY begin_time ASC',tup))
    
        #get all records from the query in python list format. 
        list=r.fetchall()

        #times are in unix time - want to use datetime instead
        output=[]

        for l in list:
            insertion_time=datetime.datetime.utcfromtimestamp(l[0])
            start_time=datetime.datetime.utcfromtimestamp(l[1])
            ref_time=datetime.datetime.utcfromtimestamp(l[2])
            end_time=datetime.datetime.utcfromtimestamp(l[3])
            event_type=l[4]
            #get a string descriptor of the event type
            event_type_info=self._lytaf_event2string(l[4])

            #create output tuple for each entry in list
            entry=(insertion_time,start_time,ref_time,end_time,event_type,event_type_info[0])
            output.append(entry)

        self.lytaf=output
        #return output
        return

    def _lytaf_event2string(self,integers):
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

    @staticmethod
    def _get_url_for_date(date):
        """Returns a URL to the LYRA data for the specified date
        """
        dt = parse_time(date or datetime.datetime.utcnow())

        # Filename
        filename = "lyra_%s000000_lev%d_%s.fits" % (dt.strftime('%Y%m%d-'),
                                                    2, 'std')
        # URL
        base_url = "http://proba2.oma.be/lyra/data/bsd/"
        url_path = urlparse.urljoin(dt.strftime('%Y/%m/%d/'), filename)
        return urlparse.urljoin(base_url, url_path)

    @classmethod
    def _get_default_uri(cls):
        """Look for and download today's LYRA data"""
        return cls._get_url_for_date(datetime.datetime.utcnow())

    @staticmethod
    def _parse_fits(filepath):
        """Loads LYRA data from a FITS file"""
        # Open file with PyFITS
        hdulist = pyfits.open(filepath)
        fits_record = hdulist[1].data
        #secondary_header = hdulist[1].header

        # Start and end dates.  Different LYRA FITS files have
        # different tags for the date obs.
        if 'date-obs' in hdulist[0].header:
            start_str = hdulist[0].header['date-obs']
        elif 'date_obs' in hdulist[0].header:
            start_str = hdulist[0].header['date_obs']
        #end_str = hdulist[0].header['date-end']

        #start = datetime.datetime.strptime(start_str, '%Y-%m-%dT%H:%M:%S.%f')
        start = parse_time(start_str)
        #end = datetime.datetime.strptime(end_str, '%Y-%m-%dT%H:%M:%S.%f')

        # First column are times
        times = [start + datetime.timedelta(0, n) for n in fits_record.field(0)]

        # Rest of columns are the data
        table = {}

        for i, col in enumerate(fits_record.columns[1:-1]):
            table[col.name] = fits_record.field(i + 1)

        # Return the header and the data
        return hdulist[0].header, pandas.DataFrame(table, index=times)
