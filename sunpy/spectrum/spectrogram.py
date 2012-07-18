# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import datetime
import urllib2

from itertools import izip
from copy import copy, deepcopy

import numpy as np
import pyfits

from bs4 import BeautifulSoup

from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

from sunpy.time import parse_time
from sunpy.util.util import to_signed
from sunpy.spectrum.spectrum import Spectrum

# This should not be necessary, as observations do not take more than a day
# but it is used for completeness' and extendibility's sake.
SECONDS_PER_DAY = 86400

# Used for COPY_PROPERTIES
REFERENCE = 0
COPY = 1
DEEPCOPY = 2


# XXX: Probably make base class.
# XXX: t_res and f_res are probably redundant because they are in
# .shape
class Spectrogram(np.ndarray):
    # Contrary to what pylint may think, this is not an old-style class.
    # pylint: disable=E1002,W0142,R0902

    # This needs to list all attributes that need to be
    # copied to maintain the object and how to handle them.
    COPY_PROPERTIES = [
        ('time_axis', COPY),
        ('freq_axis', COPY),
        ('start', REFERENCE),
        ('end', REFERENCE),
        ('_gstart', REFERENCE),
        ('_gend', REFERENCE),
        ('t_label', REFERENCE),
        ('f_label', REFERENCE),
        ('content', REFERENCE),
        ('timedelta', REFERENCE),
    ]

    @classmethod
    def _new_with_params(cls, data, params):
        """ Implementation detail. """
        obj = cls.__new__(cls, data)
        for key, value in params.iteritems():
            setattr(obj, key, value)
        return obj
    
    def slice(self, y_range, x_range):
        """ Return new spectrogram reduced to the values passed
        as slices. """
        data = super(CallistoSpectrogram, self).__getitem__([y_range, x_range])
        params = vars(self).copy()

        soffset = 0 if x_range.start is None else x_range.start
        eoffset = self.shape[1] if x_range.stop is None else x_range.stop

        fsoffset = 0 if y_range.start is None else y_range.start
        feoffset = self.shape[0] if y_range.stop is None else y_range.stop
        
        params.update({
            'time_axis': self.time_axis[x_range.start:x_range.stop:x_range.step],
            'freq_axis': self.freq_axis[y_range.start:y_range.stop:y_range.step],
            'start': self.start + soffset * self.timedelta,
            'end': self.start + eoffset * self.timedelta,
            'f_init': self.freq_axis[fsoffset], # XXX

        })
        return self._new_with_params(data, params)

    # This accepting arbitrary arguments makes it easier to subclass this.
    def __new__(cls, data, *args, **kwargs):
        return np.asarray(data).view(cls)

    def __init__(self, data, time_axis, freq_axis, start, end,
        t_label="Time", f_label="Frequency", content=""):
        # Because of how object creation works, there is no avoiding
        # unused arguments in this case.
        self.start = start
        self.end = end

        # Starting time of the whole measurement, the time axis is relative
        # to this.
        self._gstart = self.start
        self._gend = self.end

        self.t_label = t_label
        self.f_label = f_label

        self.time_axis = time_axis
        self.freq_axis = freq_axis

    def time_formatter(self, x, pos):
        """ This returns the label for the tick of value x at
        a specified pos on the time axis. """
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
        """ This returns the label for the tick of value x at
        a specified pos on the frequency axis. """
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

    def show(self, *args, **kwargs):
        self.plot(*args, **kwargs).show()

    def plot(self, overlays=[], colorbar=True, **matplotlib_args):
        # [] as default argument is okay here because it is only read.
        # pylint: disable=W0102,R0914

        figure = plt.figure(frameon=True)
        axes = figure.add_subplot(111)
        
        params = {
            'origin': 'lower',
        }
        params.update(matplotlib_args)
        # XXX: This should not be necessary.
        im = axes.imshow(np.array(self), **params)
        
        xa = axes.get_xaxis()
        ya = axes.get_yaxis()

        # Set the tick labels to be looked up in the two axis arrays.
        # If frequencies were scaled up, we need to reverse that here.
        xa.set_major_formatter(
            FuncFormatter(self.time_formatter)
        )

        ya.set_major_locator(MaxNLocator(integer=True, steps=[1, 5, 10]))
        ya.set_major_formatter(
            FuncFormatter(self.freq_formatter)
        )
        
        axes.set_xlabel(self.t_label)
        axes.set_ylabel(self.f_label)
        figure.suptitle(self.content)
        
        for tl in xa.get_ticklabels():
            tl.set_fontsize(10)
            tl.set_rotation(30)
        figure.add_axes(axes)
        if colorbar:
            figure.colorbar(im).set_label("Intensity")

        for overlay in overlays:
            figure, axes = overlay(figure, axes)
        return figure

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key, slice(None, None, None))

        if isinstance(key[0], slice) and isinstance(key[1], slice):
            return self.slice(key[0], key[1])
        elif isinstance(key[1], slice):
            return Spectrum( # XXX: Right class
                super(Spectrogram, self).__getitem__(key),
                self.time_axis[key[1].start:key[1].stop:key[1].step]
            )
        elif isinstance(key[0], slice):
            return Spectrum(
                super(Spectrogram, self).__getitem__(key),
                self.freq_axis[key[0].start:key[0].stop:key[0].step]
            )
        
        return super(Spectrogram, self).__getitem__(key)

    def clip_freq(self, minimum=None, maximum=None):
        """ Return a new spectrogram only consisting of frequencies
        in the interval [minimum, maximum]. """
        left = 0
        if maximum is not None:
            while self.freq_axis[left] > maximum:
                left += 1

        right = len(self.freq_axis) - 1

        if minimum is not None:
            while self.freq_axis[right] < minimum:
                right -= 1

        return self[left:right, :]

    def subtract_bg(self):
        """ Perform constant background subtraction. """
        # pylint: disable=E1101,E1103
        data = self.astype(to_signed(self.dtype))
        # Subtract average value from every frequency channel.
        tmp = (data - np.average(self, 1).reshape(self.shape[0], 1))
        # Get standard deviation at every point of time.
        # Need to convert because otherwise this class's __getitem__
        # is used which assumes two-dimensionality.
        sdevs = np.asarray(np.std(tmp, 0))

        # Get indices of values with lowest standard deviation.
        cand = sorted(xrange(self.shape[0]), key=lambda y: sdevs[y])
        # Only consider the best 5 %.
        realcand = cand[:max(1, int(0.05 * len(cand)))]

        # Average the best 5 %
        bg = np.average(self[:, realcand], 1)

        return self - bg.reshape(self.shape[0], 1)

    def clip(self, minimum=None, maximum=None):
        """ Clip values to be in the interval [minimum, maximum]. Any values
        greater than the maximum will be assigned the maximum, any values
        lower than the minimum will be assigned the minimum. If either is
        left out or None, do not clip at that side of the interval. """
        # pylint: disable=E1101
        if minimum is None:
            minimum = int(self.min())

        if maximum is None:
            maximum = int(self.max())

        new = self.copy()
        new[new < minimum] = minimum
        new[new > maximum] = maximum

        return new

    def normalize(self, min_, max_, dtype_=dtype('float32')):
        data = self.astype(dtype_)
        return (
            min_ + (max_ - min_) * (data - self.min()) /
            (self.max() - self.min()
        )


class LinearTimeSpectrogram(Spectrogram):
    COPY_PROPERTIES = Spectrogram.COPY_PROPERTIES + [
        ('t_delt', REFERENCE),
        ('t_init', REFERENCE),
    ]
     
    def __init__(self, data, time_axis, freq_axis, start, end,
        t_init, t_delt, t_label="Time", f_label="Frequency",
        content=""):
        super(LinearTimeSpectrogram, self).__init__(
            data, time_axis, freq_axis, start, end, t_label, f_label,
            content
        )

        self.t_init = t_init
        self.t_delt = t_delt

        self.timedelta = datetime.timedelta(seconds=self.t_delt)

    @classmethod
    def join_many(cls, spectrograms, arr=None, nonlinear=True):
        # XXX: Only load header and load contents of files
        # on demand.

        # XXX: This currently assumes all files are sampled with
        # the same sampling rate and have the same frequency
        # channels.
        specs = sorted(spectrograms, key=lambda x: x.t_init)
        data = specs[0]
        init = data.t_init
        
        if arr is None:
            arr = np.array([], dtype=max(sp.dtype for sp in specs)
        )

        size = sum(sp.t_res for sp in specs)

        xs = []
        for elem in specs[1:]:
            x = int((elem.t_init - init) / data.t_delt)
            xs.append(x)
            diff = (data.t_res - x)

            # If we leave out undefined values, we do not want to
            # add values here if x > t_res.
            if nonlinear:
                size -= max(0, diff)
            else:
                size -= diff

            init = elem.t_init

        # The non existing element after the last one starts after
        # the last one. Needed to keep implementation below sane.
        xs.append(specs[-1].shape[1])

        # We do that here so the user can pass a memory mapped
        # array if they'd like to.
        arr.resize((data.shape[0], size))
        time_axis = np.zeros((size,))
        sx = 0

        for x, elem in izip(xs, specs):
            if x > elem.shape[1]:
                if nonlinear:
                    x = elem.shape[1]
                else:
                    # If we want to stay linear, fill up the missing
                    # pixels with placeholder zeros.
                    filler = np.zeros((data.f_res, x - elem.shape[1]))
                    filler[:] = 0
                    elem = np.concatenate([elem, filler], 1)
            arr[:, sx:sx + x] = elem[:, :x]
            time_axis[sx:sx + x] = elem.time_axis[:x] + data.t_delt * sx
            
            sx += x
        params = {
            'time_axis': time_axis,
            'freq_axis': data.freq_axis,
            'start': data.start,
            'end': specs[-1].end,
            '_gstart': data._gstart,
            't_delt': data.t_delt, # XXX
            't_init': data.t_init,
            't_label': data.t_label,
            'f_label': data.f_label,
            'content': data.content,
            'timedelta': data.timedelta,
        }
        return cls._new_with_params(arr, params)

    def time_to_x(self, time):
        """ Return x-coordinate in spectrogram that corresponds to the
        passed datetime value. """
        # This is impossible for frequencies because that mapping
        # is not injective.
        # XXX: This assumes time is linear. As we read it from the
        # axes table, it might not be.
        diff = time - self.start
        diff_s = SECONDS_PER_DAY * diff.days + diff.seconds
        td_s = SECONDS_PER_DAY * self.timedelta.days  + self.timedelta.seconds
        # time - start = k * timedelta
        k = diff_s / td_s
        return round(k * self.shape[1])
