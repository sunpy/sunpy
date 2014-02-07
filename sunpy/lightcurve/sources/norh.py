import datetime
import os
from astropy.io import fits
from sunpy.time import parse_time

__all__ = ['NoRHLightCurve']

class NoRHLightCurve(LightCurve):

def _get_url_for_date_range():
    baseurl='ftp://solar-pub.nao.ac.jp/pub/nsro/norh/data/tcx/'
    #date is a datetime object
    year=date.strftime('%Y')
    year_trim=date.strftime('%y')
    mon=date.strftime('%m')
    day=date.strftime('%d')
    final_url_17=os.path.join(baseurl,year,month,'tca'+year_trim+mon+day)
    final_url_34=os.path.join(baseurl,year,month,'tcz'+year_trim+mon+day)

    return final_url_17

    #now want to download from these urls using the sunpy net functions

def read_norh(date):
    #then find and open the files
    hdulist=fits.open('tca20110607')
    header=hdulist[0].header
    data=hdulist[0].data

    #construct the time array from the FITS header
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

    
