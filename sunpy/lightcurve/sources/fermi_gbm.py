"""Provides programs to process and analyse Fermi/GBM lightcurve data."""

from __future__ import absolute_import

import datetime
import urlparse
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
import pandas

from sunpy.lightcurve import LightCurve
from sunpy.time import parse_time
from sunpy.util.odict import OrderedDict


__all__ = ['GBMLightCurve']

class GBMLightCurve(LightCurve):
    """
    Nobeyama Radioheliograph LightCurve.

    Examples
    --------
    >>> import sunpy
    
    >>> norh = sunpy.lightcurve.NoRHLightCurve.create('~/Data/norh/tca110607')
    >>> norh = sunpy.lightcurve.NoRHLightCurve.create('2011/08/10')
    >>> norh = sunpy.lightcurve.NoRHLightCurve.create('2011/08/10',wavelength='34')
    >>> norh.peek()

    References
    ----------
    | http://solar.nro.nao.ac.jp/norh/
    """

    def peek(self, **kwargs):
        """Plots the NoRH lightcurve"""
        plt.figure()
        axes = plt.gca()
        #data_lab=self.meta['OBS-FREQ'][0:2] + ' ' + self.meta['OBS-FREQ'][2:5]
        axes.plot(self.data.index,self.data,label=data_lab)
        axes.set_yscale("log")
        axes.set_ylim(1e-4,1)
        axes.set_title('Fermi GBM')
        axes.set_xlabel('Start time: ' + self.data.index[0].strftime('%Y-%m-%d %H:%M:%S UT'))
        axes.set_ylabel('Counts')
        axes.legend()
        plt.show()

    @classmethod
    def _get_url_for_date(cls,date, **kwargs):
        """This method retrieves the url for Fermi/GBM data for the given date."""
        baseurl='http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'
        #date is a datetime object
        final_url=urlparse.urljoin(baseurl, date.strftime('%Y/%m/%d/' + 'current/' + 'glg_cspec_n1_%y%m%d_v00.pha'))
        print final_url
        
        return final_url
     
    
    @staticmethod
    def _parse_fits(filepath):
        """This method parses GBM CSPEC data files to create summary lightcurves."""
        hdulist=fits.open(filepath)
        header=OrderedDict(hdulist[0].header)
        #these GBM files have three FITS extensions.
        #extn1 - this gives the energy range for each of the 128 energy bins
        #extn2 - this contains the data, e.g. counts, exposure time, time of observation
        #extn3 - ???
        energy_bins=hdulist[1].data
        count_data=hdulist[2].data
        misc=hdulist[3].data

        #rebin the 128 energy channels into some summary ranges
        #4-15 keV, 15 - 25 keV, 25-50 keV, 50-100 keV, 100-300 keV, 300-800 keV, 800 - 2000 keV
        #put the data in the units of counts/s/keV
        summary_counts=_bin_data_for_summary(energy_bins,count_data)
       
        #times for GBM are in Mission Elapsed Time (MET). The reference time for this is 2001-Jan-01 00:00.
        met_ref_time = parse_time('2001-01-01 00:00')
        offset_from_utc = (met_ref_time - parse_time('1979-01-01 00:00')).total_seconds()

        gbm_times=[]
        #get the time information in datetime format with the correct MET adjustment
        for t in count_data['time']:
            gbm_times.append(parse_time(t + offset_from_utc))

        column_labels=['4-15 keV','15-25 keV','25-50 keV','50-100 keV','100-300 keV','300-800 keV','800-2000 keV']
        return header, pandas.DataFrame(summary_counts, columns=column_labels, index=gbm_times)


def _bin_data_for_summary(energy_bins,count_data):

    #find the indices corresponding to some standard summary energy bins
    ebands=[4,15,25,50,100,300,800,2000]
    indices=[]
    for e in ebands:
        indices.append(np.searchsorted(energy_bins['e_max'],e))
        
    #get energy indices for summary count bins
    #indices=_get_summary_bin_indices(energy_bins)
    #print indices

    #rebin the 128 energy channels into some summary ranges
    #4-15 keV, 15 - 25 keV, 25-50 keV, 50-100 keV, 100-300 keV, 300-800 keV, 800 - 2000 keV
    #put the data in the units of counts/s/keV
    summary_counts=[]
    for i in range(0,len(count_data['counts'])):
        summary_counts.append([np.sum(count_data['counts'][i][0:indices[1]]) /
                               (count_data['exposure'][i] * (energy_bins['e_max'][indices[1]] - energy_bins['e_min'][0])),
                                np.sum(count_data['counts'][i][indices[1]:indices[2]]) /
                                (count_data['exposure'][i] * (energy_bins['e_max'][indices[2]] - energy_bins['e_min'][indices[1]])),
                               np.sum(count_data['counts'][i][indices[2]:indices[3]]) /
                                 (count_data['exposure'][i] * (energy_bins['e_max'][indices[3]] - energy_bins['e_min'][indices[2]])),
                               np.sum(count_data['counts'][i][indices[3]:indices[4]]) /
                                 (count_data['exposure'][i] * (energy_bins['e_max'][indices[4]] - energy_bins['e_min'][indices[3]])),
                               np.sum(count_data['counts'][i][indices[4]:indices[5]]) /
                                (count_data['exposure'][i] * (energy_bins['e_max'][indices[5]] - energy_bins['e_min'][indices[4]])),
                               np.sum(count_data['counts'][i][indices[5]:indices[6]]) /
                                (count_data['exposure'][i] * (energy_bins['e_max'][indices[6]] - energy_bins['e_min'][indices[5]])),
                               np.sum(count_data['counts'][i][indices[6]:indices[7]]) /
                                (count_data['exposure'][i] * (energy_bins['e_max'][indices[7]] - energy_bins['e_min'][indices[6]]))])

    return summary_counts
