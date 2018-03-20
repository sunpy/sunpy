"""Provides programs to process and analyse Fermi/GBM lightcurve data."""

from __future__ import absolute_import, print_function


from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
import pandas

from sunpy.io.fits import fits
from sunpy.instr import fermi
from sunpy.lightcurve import LightCurve

from sunpy.extern.six.moves import urllib

__all__ = ['GBMSummaryLightCurve']


class GBMSummaryLightCurve(LightCurve):
    """
    Fermi/GBM Summary Lightcurve.

    The Gamma-ray Burst Monitor (GBM) is an instrument aboard Fermi. It is meant
    to detect gamma-ray bursts but also detects solar flares. It consists of
    12 Sodium Iodide (NaI) scintillation detectors and 2 Bismuth Germanate (BGO)
    scintillation detectors. The NaI detectors cover from a few keV to about 1 MeV
    and provide burst triggers and locations. The BGO detectors cover the energy range from
    about 150 keV to about 30 MeV.

    This summary lightcurve makes use of the CSPEC (daily version) data set which
    consists of the counts accumulated every 4.096 seconds in 128 energy channels
    for each of the 14 detectors. Note that the data is re-binned from the
    original 128 into the following 8 pre-determined energy channels.

        * 4-15 keV
        * 15-25 keV
        * 25-50 keV
        * 50-100 keV
        * 100-300 keV
        * 300-800 keV
        * 800-2000 keV

    Examples
    --------
    >>> from sunpy.lightcurve import GBMSummaryLightCurve
    >>> gbm = GBMSummaryLightCurve.create('2011-06-07')   # doctest: +REMOTE_DATA
    No detector specified. Detector with smallest mean angle to Sun is n5
    Using Detector n5
    For Fermi detector pointing information, use tools in sunpy/instr/fermi
    >>> gbm.peek()   # doctest: +SKIP

    References
    ----------
    * `Fermi Mission Homepage <https://fermi.gsfc.nasa.gov>`_
    * `Fermi GBM Homepage <https://fermi.gsfc.nasa.gov/science/instruments/gbm.html>`_
    * `Fermi Science Support Center <https://fermi.gsfc.nasa.gov/ssc/>`_
    * `Fermi Data Product <https://fermi.gsfc.nasa.gov/ssc/data/access/>`_
    * `GBM Instrument Papers <http://gammaray.nsstc.nasa.gov/gbm/publications/instrument_journal_gbm.html>`_
    """

    def peek(self, **kwargs):
        """Plots the GBM lightcurve. An example can be seen below.

        .. plot::

            from sunpy.lightcurve import GBMSummaryLightCurve
            from sunpy.data.sample import GBM_TIMESERIES
            gbm = GBMSummaryLightCurve.create(GBM_TIMESERIES)
            gbm.peek()

        Parameters
        ----------
        **kwargs : dict
            Any additional plot arguments that should be used
            when plotting.
        """
        figure=plt.figure()
        axes = plt.gca()
        data_lab = self.data.columns.values

        for d in data_lab:
            axes.plot(self.data.index, self.data[d], label=d)

        axes.set_yscale("log")
        axes.set_title('Fermi GBM Summary data ' + self.meta['DETNAM'])
        axes.set_xlabel('Start time: ' +
                        self.data.index[0].strftime('%Y-%m-%d %H:%M:%S UT'))
        axes.set_ylabel('Counts/s/keV')
        axes.legend()
        figure.autofmt_xdate()

        plt.show()

    @classmethod
    def _get_url_for_date(cls, date, **kwargs):
        """Returns the url for Fermi/GBM data for the given date."""
        baseurl = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/'
        # date is a datetime object
        if 'detector' in kwargs:
            det = _parse_detector(kwargs['detector'])
            final_url = urllib.parse.urljoin(
                baseurl, date.strftime('%Y/%m/%d/' + 'current/' +
                                       'glg_cspec_' + det + '_%y%m%d_v00.pha'))
        else:
            # if user doesn't specify a detector, find the one pointing
            # closest to the Sun.'
            # OR: maybe user should have to specify detector or fail.
            det = cls._get_closest_detector_for_date(date)
            print('No detector specified. Detector with smallest mean angle '
                  'to Sun is ' + str(det))
            print('Using Detector ' + str(det))
            print('For Fermi detector pointing information, use tools in '
                  'sunpy/instr/fermi')
            final_url = urllib.parse.urljoin(
                baseurl, date.strftime('%Y/%m/%d/' + 'current/' +
                                       'glg_cspec_' + det + '_%y%m%d_v00.pha'))

        return final_url

    @classmethod
    def _get_closest_detector_for_date(cls, date, **kwargs):
        """Returns the GBM detector with the smallest mean angle to the Sun
        for the given date"""
        pointing_file = fermi.download_weekly_pointing_file(date)
        det_angles = fermi.get_detector_sun_angles_for_date(date, pointing_file)
        det_angle_means = []
        for n in det_angles.keys():
            if not n == 'time':
                det_angle_values = []
                for angle in det_angles[n]:
                    det_angle_values.append(angle.value)

                det_angle_means.append(np.mean(det_angle_values))

        best_det = 'n' +str(np.argmin(det_angle_means))
        return best_det


    @staticmethod
    def _parse_fits(filepath):
        """Parses GBM CSPEC data files to create summary lightcurves."""
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

        gbm_times=[]
        #get the time information in datetime format with the correct MET adjustment
        for t in count_data['time']:
            gbm_times.append(fermi.met_to_utc(t))
        column_labels=['4-15 keV','15-25 keV','25-50 keV','50-100 keV','100-300 keV',
                       '300-800 keV','800-2000 keV']
        return header, pandas.DataFrame(summary_counts, columns=column_labels, index=gbm_times)


def _bin_data_for_summary(energy_bins,count_data):
    """Missing doc string"""
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
                                (count_data['exposure'][i] * (energy_bins['e_max'][indices[j]] -
                                                              energy_bins['e_min'][indices[j-1]])))

        summary_counts.append(counts_in_bands)

    return summary_counts


def _parse_detector(detector):
    """Missing Doc String"""
    oklist=['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','n10','n11']
    altlist = [str(i) for i in range(12)]
    if detector in oklist:
        return detector
    elif detector in altlist:
        return 'n'+detector
    else:
        raise ValueError('Detector string could not be interpreted')
