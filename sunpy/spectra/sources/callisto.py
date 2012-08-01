# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import os
import datetime
import urllib2

import numpy as np
import pyfits

from collections import defaultdict

from bs4 import BeautifulSoup

from sunpy.time import parse_time
from sunpy.util.magicfunc import MagicFunc
from sunpy.spectra.spectrogram import LinearTimeSpectrogram, REFERENCE

TIME_STR = "%Y%m%d%H%M%S"
DEFAULT_URL = 'http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto/'
_DAY = datetime.timedelta(days=1)

def _buffered_write(inp, outp, buffer_size):
    """ Implementation detail. Write from inp to outp in chunks of
    buffer_size. """
    while True:
        read = inp.read(buffer_size)
        if not read:
            break
        outp.write(read)


def query(start, end, instruments=None, url=DEFAULT_URL):
    """ Get URLs for callisto data from instruments between start and end.
    
    Parameters
    ----------
    start : parse_time compatible
    end : parse_time compatible
    instruments : sequence
        Sequence of instruments whose data is requested.
    url : str
        Base URL for the request.
    """
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
                # If the split fails, the file name does not match out format,
                # so we skip it and continue to the next iteration of the loop.
                continue
            point = datetime.datetime.strptime(date + time, TIME_STR)
            opn.close()
            if (instruments is not None and
                inst not in instruments and 
                (inst, int(no)) not in instruments):
                continue
            
            if start <= point <= end:
                yield directory + href
        day += _DAY


def download(urls, directory):
    """ Download files from urls into directory.
    
    Parameters
    ----------
    urls : list of str
        urls of the files to retrieve
    directory : str
        directory to save them in
    """
    paths = []
    for url in urls:
        _, filename = os.path.split(url)
        path = os.path.join(directory, filename)
        fd = open(path, 'w')
        src = urllib2.urlopen(url)
        _buffered_write(src, fd, 4096)
        fd.close()
        src.close()
        paths.append(path)
    return paths


def parse_header_time(date, time):
    """ Return datetime object from date and time fields of header. """
    if time is not None:
        date = date + 'T' + time
    return parse_time(date)


class CallistoSpectrogram(LinearTimeSpectrogram):
    """ Classed used for dynamic spectra coming from the Callisto network.
    
    
    Additional (not inherited) parameters
    -------------------------------------
    header : pyfits.Header
        main header of the FITS file
    axes_header : pyfits.Header
        header foe the axes table
    swapped : boolean
        flag that specifies whether originally in the file the x-axis was
        frequency
    """
    make = MagicFunc()
    # Contrary to what pylint may think, this is not an old-style class.
    # pylint: disable=E1002,W0142,R0902

    # This needs to list all attributes that need to be
    # copied to maintain the object and how to handle them.
    COPY_PROPERTIES = LinearTimeSpectrogram.COPY_PROPERTIES + [
        ('header', REFERENCE),
        ('swapped', REFERENCE),
        ('axes_header', REFERENCE)
    ]
    
    # List of instruments retrieved in July 2012 from
    # http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto/
    INSTRUMENTS = set([
        'ALASKA', 'ALMATY', 'BIR', 'DARO', 'HB9SCT', 'HUMAIN',
        'HURBANOVO', 'KASI', 'KENYA', 'KRIM', 'MALAYSIA', 'MRT1',
        'MRT2', 'OOTY', 'OSRA', 'SWMC', 'TRIEST', 'UNAM'
    ])

    def save(self, filepath):
        """ Save modified spectrogram back to filepath.
        
        Parameters
        ----------
        filepath : str
            path to save the spectrogram to
        """
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

    @classmethod
    def read(cls, filename, **kwargs):
        """ Read in FITS file and return a new CallistoSpectrogram. 
        Any unknown (i.e. any except filename) keyword arguments get
        passed to pyfits.open.
        
        Parameters
        ----------
        filename : str
            path of the file to read
        """
        fl = pyfits.open(filename, **kwargs)
        data = fl[0].data
        axes = fl[1]
        header = fl[0].header

        start = parse_header_time(
            header['DATE-OBS'], header.get('TIME-OBS', header.get('TIME$_OBS'))
        )
        end = parse_header_time(
            header['DATE-END'], header.get('TIME-END', header.get('TIME$_END'))
        )

        swapped = "time" not in header["CTYPE1"].lower()
        
        # Swap dimensions so x-axis is always time.
        if swapped:
            t_delt = header["CDELT2"]
            t_init = header["CRVAL2"] - t_delt * header["CRPIX2"]
            t_label = header["CTYPE2"]

            f_delt = header["CDELT1"]
            f_init = header["CRVAL1"] - t_delt * header["CRPIX1"]
            f_label = header["CTYPE1"]
            data = data.transpose()
        else:
            t_delt = header["CDELT1"]
            t_init = header["CRVAL1"] - t_delt * header["CRPIX1"]
            t_label = header["CTYPE1"]

            f_delt = header["CDELT2"]
            f_init = header["CRVAL2"] - t_delt * header["CRPIX2"]
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
                np.linspace(0, data.shape[1] - 1) * t_delt + t_init # pylint: disable=E1101

        if fq is not None:  
            freq_axis = np.squeeze(fq)
        else:
            freq_axis = \
                np.linspace(0, data.shape[0] - 1) * f_delt + f_init # pylint: disable=E1101

        content = header["CONTENT"].split(" ", 1)[1]
        content = start.strftime("%d %b %Y")+ ' ' + content
        return cls(data, time_axis, freq_axis, start, end, t_init, t_delt,
            t_label, f_label, content, header, axes.header, swapped)

        
    def __init__(self, data, time_axis, freq_axis, start, end,
            t_init, t_delt, t_label, f_label, content, header, axes_header,
            swapped):
        # Because of how object creation works, there is no avoiding
        # unused arguments in this case.
        # pylint: disable=W0613
        
        super(CallistoSpectrogram, self).__init__(
            data, time_axis, freq_axis, start, end,
            t_init, t_delt, t_label, f_label,
            content
        )

        self.header = header
        self.axes_header = axes_header
        self.swapped = swapped

    @classmethod
    def is_datasource_for(cls, header):
        """ Check if class supports data from the given FITS file.
        
        Parameters
        ----------
        header : pyfits.Header
            main header of the FITS file
        """
        return header.get('instrument', '').strip() in cls.INSTRUMENTS

    def remove_border(self):
        """ Remove duplicate entries on the borders. """
        left = 0
        while self.freq_axis[left] == self.freq_axis[0]:
            left += 1
        right = self.shape[0] - 1
        while self.freq_axis[right] == self.freq_axis[-1]:
            right -= 1
        return self[left-1:right+2, :]

    @classmethod
    def read_many(cls, filenames, sort_by=None):
        """ Return list of CallistoSpectrogram objects read from filenames.
        
        Parameters
        ----------
        filenames : list of str
            list of paths to read from
        sort_by : str
            optional attribute of the resulting objects to sort from, e.g.
            start to sort by starting time.
        """
        objs = map(cls.read, filenames)
        if sort_by is not None:
            objs.sort(key=lambda x: getattr(x, sort_by))
        return objs

    @classmethod
    def from_file(cls, filename):
        return cls.read(filename)
    
    @classmethod
    def from_dir(cls, directory):
        return cls.read_many(
            os.path.join(directory, elem) for elem in os.listdir(directory)
        )
    
    @classmethod
    def from_url(cls, url):
        """ Return CallistoSpectrogram read from URL.
        
        Parameters
        ----------
        url : str
            URL to retrieve the data from
        """
        return cls.read(url)
    
    @classmethod
    def from_range(cls, instrument, start, end):
        """ Automatically download data from instrument between start and
        end and join it together.
        
        Parameters
        ----------
        instrument : str
            instrument to retrieve the data from
        start : parse_time compatible
            start of the measurement
        end : parse_time compatible
            end of the measurement
        """
        start = parse_time(start)
        end = parse_time(end)
        urls = query(start, end, [instrument])
        data = map(cls.from_url, urls)
        freq_buckets = defaultdict(list)
        for elem in data:
            freq_buckets[tuple(elem.freq_axis)].append(elem)
        return cls.combine_frequencies(
            [cls.join_many(elem) for elem in freq_buckets.itervalues()]
        )


CallistoSpectrogram.make.add(
    CallistoSpectrogram.from_file,
    lambda filename: os.path.isfile(filename)
)
CallistoSpectrogram.make.add(
    CallistoSpectrogram.from_dir,
    lambda directory: os.path.isdir(directory)
)
CallistoSpectrogram.make.add(
    CallistoSpectrogram.from_url
)
CallistoSpectrogram.make.add(
    CallistoSpectrogram.from_range
)


if __name__ == "__main__":
    opn = CallistoSpectrogram.read("callisto/BIR_20110922_103000_01.fit")
    opn.subtract_bg().clip(0).plot(ratio=2).show()
    print "Press return to exit"
    raw_input()
