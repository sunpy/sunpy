"""Provides programs to process and analyse Fermi/GBM lightcurve data."""

from __future__ import absolute_import

import datetime
import urlparse
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import pandas
import warnings

from sunpy.instr import fermi
from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time
from sunpy.util.odict import OrderedDict


__all__ = ['GBMSummaryLightCurve']

class GBMSummaryLightCurve(LightCurve):
    """
    Fermi/GBM Summary Lightcurve.

    Examples
    --------
    >>> import sunpy
    
    >>> gbm = sunpy.lightcurve.GBMLightCurve.create('2011-06-07')
    >>> gbm.peek()

    References
    ----------
    | http://gammaray.nsstc.nasa.gov/gbm/
    """

    def peek(self, **kwargs):
        """Plots the GBM lightcurve"""
        figure=plt.figure()
        axes = plt.gca()
        data_lab=self.data.columns.values

        for d in data_lab:
            axes.plot(self.data.index,self.data[d],label=d)
        
        axes.set_yscale("log")
        axes.set_title('Fermi GBM Summary data ' + self.meta['DETNAM'])
        axes.set_xlabel('Start time: ' + self.data.index[0].strftime('%Y-%m-%d %H:%M:%S UT'))
        axes.set_ylabel('Counts/s/keV')
        axes.legend()
        figure.autofmt_xdate()
       
        plt.show()

    @classmethod
    def _get_url_for_date(cls,date, **kwargs):
        """This method retrieves the url for Fermi/GBM data for the given date."""
        baseurl='http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'
        #date is a datetime object
        if 'detector' in kwargs:
            det=_parse_detector(kwargs['detector'])
            final_url=urlparse.urljoin(baseurl, date.strftime('%Y/%m/%d/' + 'current/' + 'glg_cspec_'+det+'_%y%m%d_v00.pha'))
        else:
            #if user doesn't specify a detector, find the one pointing closest to the Sun.'
            #OR: maybe user should have to specify detector or fail.
            det = cls._get_closest_detector_for_date(date)
            print 'No detector specified. Detector with smallest mean angle to Sun is ' + str(det)
            print 'Using Detector ' + str(det)
            print 'For Fermi detector pointing information, use tools in sunpy/instr/fermi' 
            final_url=urlparse.urljoin(baseurl, date.strftime('%Y/%m/%d/' + 'current/' + 'glg_cspec_'+det+'_%y%m%d_v00.pha'))
        
        return final_url

    @classmethod
    def _get_closest_detector_for_date(cls,date,**kwargs):
        '''This method returns the GBM detector with the smallest mean angle to the Sun for the given date'''
        pointing_file = fermi.download_weekly_pointing_file(date)
        det_angles = fermi.get_detector_sun_angles_for_date(date,pointing_file)
        det_angle_means=[]
        for n in det_angles.keys():
            if not n == 'time':
                det_angle_values=[]
                for angle in det_angles[n]:
                    det_angle_values.append(angle.value)
                    
                det_angle_means.append(np.mean(det_angle_values)) 
        
       # for d in det_angles:
        #    det_angle_means.append(np.mean(d))
        best_det = 'n' +str(np.argmin(det_angle_means))
        return best_det
     
    
    @staticmethod
    def _parse_fits(filepath):
        """This method parses GBM CSPEC data files to create summary lightcurves."""
        hdulist=fits.open(filepath)
        header=OrderedDict(hdulist[0].header)
        #these GBM files have three FITS extensions.
        #extn1 - this gives the energy range for each of the 128 energy bins
        #extn2 - this contains the data, e.g. counts, exposure time, time of observation
        #extn3 - eclipse times?
        energy_bins=hdulist[1].data
        count_data=hdulist[2].data
        misc=hdulist[3].data


        #rebin the 128 energy channels into some summary ranges
        #4-15 keV, 15 - 25 keV, 25-50 keV, 50-100 keV, 100-300 keV, 300-800 keV, 800 - 2000 keV
        #put the data in the units of counts/s/keV
        summary_counts=_bin_data_for_summary(energy_bins,count_data)
       
        #times for GBM are in Mission Elapsed Time (MET). The reference time for this is 2001-Jan-01 00:00.
        #met_ref_time = parse_time('2001-01-01 00:00')
        #offset_from_utc = (met_ref_time - parse_time('1979-01-01 00:00')).total_seconds()

        gbm_times=[]
        #get the time information in datetime format with the correct MET adjustment
        for t in count_data['time']:
            #gbm_times.append(parse_time(t + offset_from_utc))
            gbm_times.append(fermi.met_to_utc(t))
        column_labels=['4-15 keV','15-25 keV','25-50 keV','50-100 keV','100-300 keV','300-800 keV','800-2000 keV']
        return header, pandas.DataFrame(summary_counts, columns=column_labels, index=gbm_times)


def _bin_data_for_summary(energy_bins,count_data):

    #find the indices corresponding to some standard summary energy bins
    ebands=[4,15,25,50,100,300,800,2000]
    indices=[]
    for e in ebands:
        indices.append(np.searchsorted(energy_bins['e_max'],e))
        
    #rebin the 128 energy channels into some summary ranges
    #4-15 keV, 15 - 25 keV, 25-50 keV, 50-100 keV, 100-300 keV, 300-800 keV, 800 - 2000 keV
    #put the data in the units of counts/s/keV
    summary_counts=[]
    for i in range(0,len(count_data['counts'])):
        counts_in_bands=[]
        for j in range(1,len(ebands)):
            counts_in_bands.append(np.sum(count_data['counts'][i][indices[j-1]:indices[j]]) /
                                (count_data['exposure'][i] * (energy_bins['e_max'][indices[j]] - energy_bins['e_min'][indices[j-1]])))
            
        summary_counts.append(counts_in_bands)

    return summary_counts

   
def _parse_detector(detector):
    oklist=['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','n10','n11']
    altlist = [str(i) for i in range(12)]
    if detector in oklist:
        return detector
    elif detector in altlist:
        return 'n'+detector
    else:
        raise ValueError('Detector string could not be interpreted')
    
        
        

