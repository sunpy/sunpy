# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import datetime

from copy import copy, deepcopy

import numpy as np
import pyfits

from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

from sunpy.time import parse_time

# This should not be necessary, as observations do not take more than a day
# but it is used for completeness' and extendibility's sake.
SECONDS_PER_DAY = 86400

# Used for COPY_PROPERTIES
REFERENCE = 0
COPY = 1
DEEPCOPY = 2


def repeat_lines(data, times):
    """ Simple lossless scaling method for integer factors """
    new = np.zeros((times * data.shape[0], data.shape[1]))
    for line in xrange(data.shape[0]):
        new[times * line: times * line + times, :] = data[line, :]
    return new


def parse_header_time(date, time):
    if time is not None:
        date = date + 'T' + time
    return parse_time(date)


# XXX: Probably make base class.
class CallistoSpectrogram(np.ndarray):
    # Contrary to what pylint may think, this is not an old-style class.
    # pylint: disable=E1002,W0142,R0902

    # This needs to list all attributes that need to be
    # copied to maintain the object and how to handle them.
    COPY_PROPERTIES = [
        ('header', REFERENCE),
        ('time_axis', COPY),
        ('freq_axis', COPY),
        ('start', REFERENCE),
        ('end', REFERENCE),
        ('_gstart', REFERENCE),
        ('_gend', REFERENCE),
        ('t_delt', REFERENCE),
        ('t_init', REFERENCE),
        ('t_label', REFERENCE),
        ('t_res', REFERENCE),
        ('f_delt', REFERENCE),
        ('f_init', REFERENCE),
        ('f_label', REFERENCE),
        ('f_res', REFERENCE),
        ('content', REFERENCE)
    ]

    def save(self, filepath):
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
        header = self.header.copy()

        if self.swapped:
            header['NAXIS2'] = self.t_res
            header['NAXIS1'] = self.f_res
        else:
            header['NAXIS1'] = self.t_res
            header['NAXIS2'] = self.f_res
        return header

    def __new__(cls, data, axes=None, header=None):
        # pylint: disable=W0613
        if header is not None:
            # Always put time on the x-axis.
            if "time" not in header["CTYPE1"].lower():
                data = data.transpose()
        
        return np.asarray(data).view(cls)
    
    @classmethod
    def _new_with_params(cls, data, params):
        obj = cls.__new__(cls, data)
        for key, value in params.iteritems():
            setattr(obj, key, value)
        return obj
    
    @classmethod
    def read(cls, filename):
        fl = pyfits.open(filename)
        return cls(fl[0].data, fl[1], fl[0].header)
    
    def slice(self, y_range, x_range):
        params = slice(y_range[0], y_range[1]), slice(x_range[0], x_range[1])
        data = super(CallistoSpectrogram, self).__getitem__(params)
        params = vars(self).copy()

        soffset = 0 if x_range[0] is None else x_range[0]
        eoffset = self.t_res if x_range[1] is None else x_range[1]

        fsoffset = 0 if y_range[0] is None else y_range[0]
        feoffset = self.f_res if y_range[1] is None else y_range[1]
        
        params.update({
            'time_axis': self.time_axis[x_range[0]:x_range[1]],
            'freq_axis': self.freq_axis[y_range[0]:y_range[1]],
            'start': self.start + soffset * self.timedelta,
            'end': self.start + eoffset * self.timedelta,
            'f_init': self.freq_axis[fsoffset],
            'f_res': feoffset - fsoffset,
            't_res': eoffset - soffset,
        })
        return self._new_with_params(data, params)
        
    def __init__(self, data, axes, header):
        # Because of how object creation works, there is no avoiding
        # unused arguments in this case.
        # pylint: disable=W0613

        self.header = header
        self.content = header["CONTENT"]

        self.axes_header = axes.header
        self.axes = axes
        
        self.start = parse_header_time(
            header['DATE-OBS'], header.get('TIME-OBS', header.get('TIME$_OBS'))
        )
        self.end = parse_header_time(
            header['DATE-END'], header.get('TIME-END', header.get('TIME$_END'))
        )

        # Starting time of the whole measurement, the time axis is relative
        # to this.
        self._gstart = self.start
        self._gend = self.end

        self.swapped = "time" not in header["CTYPE1"].lower()
        
        if self.swapped:
            self.t_res = header["NAXIS2"]
            self.t_delt = header["CDELT2"]
            self.t_init = header["CRVAL2"]
            self.t_label = header["CTYPE2"]

            self.f_res = header["NAXIS1"]
            self.f_delt = header["CDELT1"]
            self.f_init = header["CRVAL1"]
            self.f_label = header["CTYPE1"]  
        else:
            self.t_res = header["NAXIS1"]
            self.t_delt = header["CDELT1"]
            self.t_init = header["CRVAL1"]
            self.t_label = header["CTYPE1"]

            self.f_res = header["NAXIS2"]
            self.f_delt = header["CDELT2"]
            self.f_init = header["CRVAL2"]
            self.f_label = header["CTYPE2"]

        diff = self.end - self.start
        d_secs = diff.days * SECONDS_PER_DAY + diff.seconds

        # In principle CDELT1, but unrounded.        
        self.t_delt = d_secs / float(self.t_res - 1)
        self.timedelta = datetime.timedelta(
            seconds=d_secs / float(self.t_res - 1)
        )
        
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
            self.time_axis = np.squeeze(tm)
        else:
            # Otherwise, assume it's linear.
            self.time_axis = \
                np.linspace(0, self.t_res - 1) * self.t_delt + self.t_init

        if fq is not None:  
            self.freq_axis = np.squeeze(fq)
        else:   
            self.freq_axis = \
                np.linspace(0, self.f_res - 1) * self.f_delt + self.f_init

    def time_formatter(self, x, pos):
        """ This returns the label for the tick of value x at
        a specified pos on the axis. """
        # Callback, cannot avoid unused arguments.
        # pylint: disable=W0613
        try:
            return self.format_time(
                self._gstart + datetime.timedelta(
                    seconds=self.time_axis[int(x)]
                )
            )
        except IndexError:
            return None

    def freq_formatter(self, x, pos):
        # Callback, cannot avoid unused arguments.
        # pylint: disable=W0613
        try:
            return self.format_freq(self.freq_axis[x])
        except IndexError:
            return None 

    def __array_finalize__(self, obj):
        if self is obj:
            return

        for prop, cpy in self.COPY_PROPERTIES:
            elem = getattr(obj, prop, None)
            if cpy == 1:
                elem = copy(elem)
            if cpy == 2:
                elem = deepcopy(elem)

            setattr(self, prop, elem)
    
    @staticmethod
    def format_time(time):
        """ Override to configure default plotting """
        return time.strftime("%H:%M:%S")
    
    @staticmethod
    def format_freq(freq):
        """ Override to configure default plotting """
        return "%.2f" % freq

    def good_ratio(self, ratio):
        # pylint: disable=E1101
        return self.shape[1] % (ratio * self.shape[0]) == 0

    def plot(self, overlays=[], colorbar=True, ratio=None, **matplotlib_args):
        # [] as default argument is okay here because it is only read.
        # pylint: disable=W0102
        if ratio is not None:
            size = self.shape[1] / ratio # pylint: disable=E1101
            times = size / self.shape[0] # pylint: disable=E1101
            data = repeat_lines(self, times)
        else:
            times = 1
            data = np.array(self)

        figure = plt.figure(frameon=True)
        axes = figure.add_subplot(111)
        
        params = {
            'origin': 'lower',
        }
        params.update(matplotlib_args)
        im = axes.imshow(data, **params)
        
        xa = axes.get_xaxis()
        ya = axes.get_yaxis()

        xa.set_major_formatter(
            FuncFormatter(lambda x, pos: self.time_formatter(x // times, pos))
        )

        ya.set_major_locator(MaxNLocator(integer=True, steps=[1, 5, 10]))
        ya.set_major_formatter(
            FuncFormatter(lambda x, pos: self.freq_formatter(x // times, pos))
        )
        
        axes.set_xlabel(self.t_label)
        axes.set_ylabel(self.f_label)
        figure.suptitle(self.content)
        
        for tl in xa.get_ticklabels():
            tl.set_fontsize(10)
            tl.set_rotation(30)
        figure.add_axes(axes)
        if colorbar:
            figure.colorbar(im)

        for overlay in overlays:
            figure, axes = overlay(figure, axes)
        return figure

    def __getitem__(self, key):
        # XXX: Fix
        if isinstance(key, tuple):
            if isinstance(key[0], slice) and isinstance(key[1], slice):
                x_range = [key[1].start, key[1].stop]
                y_range = [key[0].start, key[0].stop]

                return self.slice(y_range, x_range)
        
        return super(CallistoSpectrogram, self).__getitem__(key)

    def time_to_x(self, time):
        # This is impossible for frequencies because that mapping
        # is not injective.
        # XXX: This assumes time is linear. As we read it from the
        # axes table, it might not be.
        diff = time - self.start
        diff_s = SECONDS_PER_DAY * diff.days + diff.seconds
        td_s = SECONDS_PER_DAY * self.timedelta.days  + self.timedelta.seconds
        # time - start = k * timedelta
        k = diff_s / td_s
        return round(k * self.t_res)

    def clip_freq(self, minimum=None, maximum=None):
        left = 0
        if maximum is not None:
            while self.freq_axis[left] > maximum:
                left += 1

        right = len(self.freq_axis) - 1

        if minimum is not None:
            while self.freq_axis[right] < minimum:
                right -= 1

        return self[left:right,:]

    def subtract_bg(self):
        # pylint: disable=E1101,E1103
        tmp = (self - np.average(self, 1).reshape(self.shape[0], 1))
        sdevs = np.std(tmp, 0)

        cand = sorted(xrange(self.shape[0]), key=lambda y: sdevs[y])
        realcand = cand[:max(1, int(0.05 * len(cand)))]

        bg = np.average(self[:, realcand], 1)

        return self - bg.reshape(self.shape[0], 1)

    @staticmethod
    def is_datasource_for(header):
        return header.get('instrument', '').startswith('BIR')

    def clip(self, minimum=None, maximum=None):
        # pylint: disable=E1101
        if minimum is None:
            minimum = int(self.min())

        if maximum is None:
            maximum = int(self.max())

        new = self.copy()
        new[new < minimum] = minimum
        new[new > maximum] = maximum

        return new


if __name__ == "__main__":
    opn = CallistoSpectrogram.read("callisto/BIR_20110922_103000_01.fit")
    opn.subtract_bg().clip(0).plot().show()
    print "Press return to exit"
    raw_input()
