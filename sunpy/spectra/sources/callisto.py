# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import datetime
import urllib2

import numpy as np

from astropy.io import fits

from collections import defaultdict

from bs4 import BeautifulSoup

from scipy.optimize import leastsq
from scipy.ndimage import gaussian_filter1d

from sunpy.time import parse_time
from sunpy.util import polyfun_at, minimal_pairs
from sunpy.util.cond_dispatch import ConditionalDispatch, run_cls
from sunpy.util.net import download_file

from sunpy.spectra.spectrogram import LinearTimeSpectrogram, REFERENCE


__all__ = ['CallistoSpectrogram']

TIME_STR = "%Y%m%d%H%M%S"
DEFAULT_URL = 'http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto/'
_DAY = datetime.timedelta(days=1)

DATA_SIZE = datetime.timedelta(seconds=15*60)

def parse_filename(href):
    name = href.split('.')[0]
    try:
        inst, date, time, no = name.rsplit('_')
        dstart = datetime.datetime.strptime(date + time, TIME_STR)
    except ValueError:
        # If the split fails, the file name does not match out
        # format,so we skip it and continue to the next
        # iteration of the loop.
        return None
    return inst, no, dstart


PARSERS = [
    # Everything starts with ""
    ("", parse_filename)
]
def query(start, end, instruments=None, url=DEFAULT_URL):
    """Get URLs for callisto data from instruments between start and end.

    Parameters
    ----------
    start : `~sunpy.time.parse_time` compatible
    end : `~sunpy.time.parse_time` compatible
    instruments : sequence
        Sequence of instruments whose data is requested.
    url : str
        Base URL for the request.
    """
    day = datetime.datetime(start.year, start.month, start.day)
    while day <= end:
        directory = url + day.strftime('%Y/%m/%d/')
        opn = urllib2.urlopen(directory)
        try:
            soup = BeautifulSoup(opn)
            for link in soup.find_all("a"):
                href = link.get("href")
                for prefix, parser in PARSERS:
                    if href.startswith(prefix):
                        break

                result = parser(href)
                if result is None:
                    continue
                inst, no, dstart = result

                if (instruments is not None and
                    inst not in instruments and
                    (inst, int(no)) not in instruments):
                    continue

                dend = dstart + DATA_SIZE
                if dend > start and dstart < end:
                    yield directory + href
        finally:
            opn.close()
        day += _DAY


def download(urls, directory):
    """Download files from urls into directory.

    Parameters
    ----------
    urls : list of str
        urls of the files to retrieve
    directory : str
        directory to save them in
    """
    return [download_file(url, directory) for url in urls]


def _parse_header_time(date, time):
    """Returns `~datetime.datetime` object from date and time fields of header. """
    if time is not None:
        date = date + 'T' + time
    return parse_time(date)


class CallistoSpectrogram(LinearTimeSpectrogram):
    """ Classed used for dynamic spectra coming from the Callisto network.

    Attributes
    ----------
    header : `~astropy.io.fits.Header`
        main header of the FITS file
    axes_header : `~astropy.io.fits.Header`
        header foe the axes table
    swapped : bool
        flag that specifies whether originally in the file the x-axis was
        frequency
    """
    # XXX: Determine those from the data.
    SIGMA_SUM = 75
    SIGMA_DELTA_SUM = 20
    _create = ConditionalDispatch.from_existing(LinearTimeSpectrogram._create)
    create = classmethod(_create.wrapper())
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
        data = fits.PrimaryHDU(self, header=main_header)
        ## XXX: Update axes header.

        freq_col = fits.Column(
            name="frequency", format="D8.3", array=self.freq_axis
        )
        time_col = fits.Column(
            name="time", format="D8.3", array=self.time_axis
        )
        cols = fits.ColDefs([freq_col, time_col])
        table = fits.new_table(cols, header=self.axes_header)

        hdulist = fits.HDUList([data, table])
        hdulist.writeto(filepath)

    def get_header(self):
        """Returns the updated header."""
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
        """Reads in FITS file and return a new CallistoSpectrogram.
        Any unknown (i.e. any except filename) keyword arguments get
        passed to fits.open.

        Parameters
        ----------
        filename : str
            path of the file to read
        """
        fl = fits.open(filename, **kwargs)
        data = fl[0].data
        axes = fl[1]
        header = fl[0].header

        start = _parse_header_time(
            header['DATE-OBS'], header.get('TIME-OBS', header.get('TIME$_OBS'))
        )
        end = _parse_header_time(
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

        content = header["CONTENT"]
        instruments = set([header["INSTRUME"]])

        return cls(
            data, time_axis, freq_axis, start, end, t_init, t_delt,
            t_label, f_label, content, instruments,
            header, axes.header, swapped
        )


    def __init__(self, data, time_axis, freq_axis, start, end,
            t_init=None, t_delt=None, t_label="Time", f_label="Frequency",
            content="", instruments=None, header=None, axes_header=None,
            swapped=False):
        # Because of how object creation works, there is no avoiding
        # unused arguments in this case.
        # pylint: disable=W0613

        super(CallistoSpectrogram, self).__init__(
            data, time_axis, freq_axis, start, end,
            t_init, t_delt, t_label, f_label,
            content, instruments
        )

        self.header = header
        self.axes_header = axes_header
        self.swapped = swapped

    @classmethod
    def is_datasource_for(cls, header):
        """Check if class supports data from the given FITS file.

        Parameters
        ----------
        header : `~astropy.io.fits.Header`
            main header of the FITS file
        """
        return header.get('instrume', '').strip() in cls.INSTRUMENTS

    def remove_border(self):
        """Remove duplicate entries on the borders."""
        left = 0
        while self.freq_axis[left] == self.freq_axis[0]:
            left += 1
        right = self.shape[0] - 1
        while self.freq_axis[right] == self.freq_axis[-1]:
            right -= 1
        return self[left-1:right+2, :]

    @classmethod
    def read_many(cls, filenames, sort_by=None):
        """Returns a list of CallistoSpectrogram objects read from filenames.

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
    def from_range(cls, instrument, start, end, **kwargs):
        """Automatically download data from instrument between start and
        end and join it together.

        Parameters
        ----------
        instrument : str
            instrument to retrieve the data from
        start : `~sunpy.time.parse_time` compatible
            start of the measurement
        end : `~sunpy.time.parse_time` compatible
            end of the measurement
        """
        kw = {
            'maxgap': None,
            'fill': cls.JOIN_REPEAT,
        }

        kw.update(kwargs)
        start = parse_time(start)
        end = parse_time(end)
        urls = query(start, end, [instrument])
        data = map(cls.from_url, urls)
        freq_buckets = defaultdict(list)
        for elem in data:
            freq_buckets[tuple(elem.freq_axis)].append(elem)
        try:
            return cls.combine_frequencies(
                [cls.join_many(elem, **kw) for elem in freq_buckets.itervalues()]
            )
        except ValueError:
            raise ValueError("No data found.")

    def _overlap(self, other):
        """ Find frequency and time overlap of two spectrograms. """
        one, two = self.intersect_time([self, other])
        ovl = one.freq_overlap(two)
        return one.clip_freq(*ovl), two.clip_freq(*ovl)

    @staticmethod
    def _to_minimize(a, b):
        """Function to be minimized for matching to frequency channels."""
        def _fun(p):
            if p[0] <= 0.2 or abs(p[1]) >= a.max():
                return float("inf")
            return a - (p[0] * b + p[1])
        return _fun

    def _homogenize_params(self, other, maxdiff=1):
        """
        Return triple with a tuple of indices (in self and other, respectively),
        factors and constants at these frequencies.

        Parameters
        ----------
        other : `sunpy.spectra.CallistoSpectrogram`
            Spectrogram to be homogenized with the current one.
        maxdiff : float
            Threshold for which frequencies are considered equal.
        """

        pairs_indices = [
            (x, y) for x, y, d in minimal_pairs(self.freq_axis, other.freq_axis)
            if d <= maxdiff
        ]

        pairs_data = [
            (self[n_one, :], other[n_two, :]) for n_one, n_two in pairs_indices
        ]

        # XXX: Maybe unnecessary.
        pairs_data_gaussian = [
            (gaussian_filter1d(a, 15), gaussian_filter1d(b, 15))
            for a, b in pairs_data
        ]

        # If we used integer arithmetic, we would accept more invalid
        # values.
        pairs_data_gaussian64 = np.float64(pairs_data_gaussian)
        least = [
            leastsq(self._to_minimize(a,b), [1, 0])[0]
            for a, b in pairs_data_gaussian64
        ]

        factors = [x for x, y in least]
        constants = [y for x, y in least]

        return pairs_indices, factors, constants

    def homogenize(self, other, maxdiff=1):
        """ Return overlapping part of self and other as (self, other) tuple.
        Homogenize intensities so that the images can be used with
        combine_frequencies. Note that this works best when most of the
        picture is signal, so use :py:meth:`in_interval` to select the subset
        of your image before applying this method.

        Parameters
        ----------
        other : `sunpy.spectra.CallistoSpectrogram`
            Spectrogram to be homogenized with the current one.
        maxdiff : float
            Threshold for which frequencies are considered equal.
        """
        one, two = self._overlap(other)
        pairs_indices, factors, constants = one._homogenize_params(
            two, maxdiff
        )
        # XXX: Maybe (xd.freq_axis[x] + yd.freq_axis[y]) / 2.
        pairs_freqs = [one.freq_axis[x] for x, y in pairs_indices]

        # XXX: Extrapolation does not work this way.
        # XXX: Improve.
        f1 = np.polyfit(pairs_freqs, factors, 3)
        f2 = np.polyfit(pairs_freqs, constants, 3)

        return (
            one,
            two * polyfun_at(f1, two.freq_axis)[:, np.newaxis] +
                polyfun_at(f2, two.freq_axis)[:, np.newaxis]
        )

    def extend(self, minutes=15, **kwargs):
        """Requests subsequent files from the server. If minutes is negative,
        retrieve preceding files. """
        if len(self.instruments) != 1:
            raise ValueError

        instrument = iter(self.instruments).next()
        if minutes > 0:
            data = CallistoSpectrogram.from_range(
                instrument,
                self.end, self.end + datetime.timedelta(minutes=minutes)
            )
        else:
            data = CallistoSpectrogram.from_range(
                instrument,
                self.start - datetime.timedelta(minutes=-minutes), self.start
            )

        data = data.clip_freq(self.freq_axis[-1], self.freq_axis[0])
        return CallistoSpectrogram.join_many([self, data], **kwargs)

    @classmethod
    def from_url(cls, url):
        """Returns CallistoSpectrogram read from URL.

        Parameters
        ----------
        url : str
            URL to retrieve the data from

        Returns
        -------
        newSpectrogram : `sunpy.spectra.CallistoSpectrogram`
        """
        return cls.read(url)


CallistoSpectrogram._create.add(
    run_cls('from_range'),
    lambda cls, instrument, start, end: True,
    check=False
)

CallistoSpectrogram.create.im_func.__doc__ = (
    """Create CallistoSpectrogram from given input dispatching to the
    appropriate from_* function.

Possible signatures:

""" + CallistoSpectrogram._create.generate_docs()
)

if __name__ == "__main__":
    opn = CallistoSpectrogram.read("callisto/BIR_20110922_103000_01.fit")
    opn.subtract_bg().clip(0).plot(ratio=2).show()
    print "Press return to exit"
    raw_input()
