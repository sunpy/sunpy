"""Provides programs to process and analyse NoRH lightcurve data."""

from __future__ import absolute_import

import datetime
import os
from astropy.io import fits
import pandas
import sunpy
from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time
import numpy as np

__all__ = ['NoRHLightCurve']

class NoRHLightCurve(LightCurve):

    @classmethod
    def _get_url_for_date(cls,date):
        """This method retrieves the url for NoRH correlation data for the given date."""
        baseurl='ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/'
        #date is a datetime object
        year=date.strftime('%Y')
        year_trim=date.strftime('%y')
        mon=date.strftime('%m')
        day=date.strftime('%d')
        final_url_17=os.path.join(baseurl,year,mon,'tca'+year_trim+mon+day)
        final_url_34=os.path.join(baseurl,year,mon,'tcz'+year_trim+mon+day)
        
        return final_url_34

    #now want to download from these urls using the sunpy net functions
    @staticmethod
    def _parse_fits(filepath):
        """This method parses NoRH tca and tcz correlation files."""
        hdulist=fits.open(filepath)
        header=hdulist[0].header
        #for these NoRH files, the time series data is recorded in the primary HDU
        data=hdulist[0].data

        #No explicit time array in FITS file, so construct the time array from the FITS header
        obs_start_time=parse_time(header['date-obs'] + 'T' + header['crval1'])
        length=len(data)
        cadence=np.float(header['cdelt1'])
        sec_array=np.linspace(0,length-1,(length/cadence))

        deltas=[]
        for s in sec_array:
            deltas.append(datetime.timedelta(0,s))

        norh_time=[]
        for d in deltas:
            norh_time.append(obs_start_time+d)

        return header, pandas.DataFrame(data,index=norh_time)

    
