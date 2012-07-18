# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import datetime
import urllib2


import numpy as np
import pyfits

from bs4 import BeautifulSoup

from sunpy.time import parse_time
from sunpy.spectrum.spectrogram import LinearTimeSpectrogram

# This should not be necessary, as observations do not take more than a day
# but it is used for completeness' and extendibility's sake.
SECONDS_PER_DAY = 86400

# Used for COPY_PROPERTIES
REFERENCE = 0
COPY = 1
DEEPCOPY = 2


TIME_STR = "%Y%m%d%H%M%S"
DEFAULT_URL = 'http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto/'
_DAY = datetime.timedelta(days=1)

def query(start, end, instrument=None, number=None, url=DEFAULT_URL):
    day = datetime.datetime(start.year, start.month, start.day)
    while day <= end:
        directory = url + '%d/%02d/%02d/' % (day.year, day.month, day.day)
        opn = urllib2.urlopen(directory)
        soup = BeautifulSoup(opn)
        for link in soup.find_all("a"):
            href = link.get("href")
            name = href.split('.')[0]
            try:
                inst, date, time, no = name.split('_')
            except ValueError:
                continue
            point = datetime.datetime.strptime(date + time, TIME_STR)
            if instrument is not None and instrument != inst:
                continue

            if number is not None and number != no:
                continue

            if start <= point <= end:
                yield directory + href
        day += _DAY


def parse_header_time(date, time):
    """ Return datetime object from date and time fields of header. """
    if time is not None:
        date = date + 'T' + time
    return parse_time(date)


class CallistoSpectrogram(LinearTimeSpectrogram):
    # Contrary to what pylint may think, this is not an old-style class.
    # pylint: disable=E1002,W0142,R0902

    # This needs to list all attributes that need to be
    # copied to maintain the object and how to handle them.
    COPY_PROPERTIES = LinearTimeSpectrogram.COPY_PROPERTIES + [
        ('header', REFERENCE),
        ('content', REFERENCE),
        ('axes_header', REFERENCE)
    ]

    INSTRUMENTS = set(['BIR'])

    def save(self, filepath):
        """ Save modified spectrogram back to filepath. """
        main_header = self.get_header()
        data = pyfits.PrimaryHDU(self, header=main_header)
        ## XXX: Update axes header.

        freq_col = pyfits.Column(
            name="frequency", format="D8.3", array=self.freq_axis
        )
        time_col = pyfits.Column(
            name="time", format="D8.3", array=self.time_axis
        )
        cols = pyfits.ColDefs([freq_col, time_col])
        table = pyfits.new_table(cols, header=self.axes_header)

        hdulist = pyfits.HDUList([data, table])
        hdulist.writeto(filepath)   

    def get_header(self):
        """ Return updated header. """
        header = self.header.copy()

        if self.swapped:
            header['NAXIS2'] = self.shape[1] # pylint: disable=E1101
            header['NAXIS1'] = self.shape[0] # pylint: disable=E1101
        else:
            header['NAXIS1'] = self.shape[1] # pylint: disable=E1101
            header['NAXIS2'] = self.shape[0] # pylint: disable=E1101
        return header

    def __new__(cls, data, axes=None, header=None):
        # pylint: disable=W0613
        if header is not None:
            # Always put time on the x-axis.
            if "time" not in header["CTYPE1"].lower():
                data = data.transpose()
        
        return np.asarray(data).view(cls)
    
    @classmethod
    def read(cls, filename):
        """ Read in FITS file and return a new CallistoSpectrogram. """
        fl = pyfits.open(filename)
        return cls(fl[0].data, fl[1], fl[0].header)
        
    def __init__(self, data, axes, header):
        # Because of how object creation works, there is no avoiding
        # unused arguments in this case.
        # pylint: disable=W0613
        start = parse_header_time(
            header['DATE-OBS'], header.get('TIME-OBS', header.get('TIME$_OBS'))
        )
        end = parse_header_time(
            header['DATE-END'], header.get('TIME-END', header.get('TIME$_END'))
        )

        self.swapped = "time" not in header["CTYPE1"].lower()
        
        # Swap dimensions so x-axis is always time.
        if self.swapped:
            t_delt = header["CDELT2"]
            t_init = header["CRVAL2"]
            t_label = header["CTYPE2"]

            f_delt = header["CDELT1"]
            f_init = header["CRVAL1"]
            f_label = header["CTYPE1"]  
        else:
            t_delt = header["CDELT1"]
            t_init = header["CRVAL1"]
            t_label = header["CTYPE1"]

            f_delt = header["CDELT2"]
            f_init = header["CRVAL2"]
            f_label = header["CTYPE2"]

        # Table may contain the axes data. If it does, the other way of doing
        # it might be very wrong.
        if axes is not None:
            try:
                # It's not my fault. Neither supports __contains__ nor .get
                tm = axes.data['time']
            except KeyError:
                tm = None
            try:
                fq = axes.data['frequency']
            except KeyError:
                fq = None
        
        if tm is not None:
            # Fix dimensions (whyever they are (1, x) in the first place)
            time_axis = np.squeeze(tm)
        else:
            # Otherwise, assume it's linear.
            time_axis = \
                np.linspace(0, self.shape[1] - 1) * t_delt + t_init # pylint: disable=E1101

        if fq is not None:  
            freq_axis = np.squeeze(fq)
        else:
            freq_axis = \
                np.linspace(0, self.shape[0] - 1) * f_delt + f_init # pylint: disable=E1101

        super(CallistoSpectrogram, self).__init__(
            data, time_axis, freq_axis, start, end,
            t_init, t_delt, t_label, f_label,
            header['CONTENT']
        )

        self.header = header
        self.axes_header = axes.header

    @classmethod
    def is_datasource_for(cls, header):
        """ Check if class supports data from the given FITS file. """
        return header.get('instrument', '').strip() in cls.INSTRUMENTS


if __name__ == "__main__":
    opn = CallistoSpectrogram.read("callisto/BIR_20110922_103000_01.fit")
    opn.subtract_bg().clip(0).plot(ratio=2).show()
    print "Press return to exit"
    raw_input()
