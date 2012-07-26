# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import datetime
import urllib2

from random import randint
from itertools import izip
from copy import copy, deepcopy
from math import floor

import numpy as np
import pyfits

from bs4 import BeautifulSoup

from numpy import ma

from scipy import ndimage
from scipy.stats.mstats import mode

from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

from sunpy.time import parse_time
from sunpy.util.util import to_signed
from sunpy.spectra.spectrum import Spectrum

# This should not be necessary, as observations do not take more than a day
# but it is used for completeness' and extendibility's sake.
# XXX: Leap second?
SECONDS_PER_DAY = 86400

# Used for COPY_PROPERTIES
REFERENCE = 0
COPY = 1
DEEPCOPY = 2

# Maybe move to util.
def get_day(dt):
    return datetime.datetime(dt.year, dt.month, dt.day)

# XXX: Find out why imshow(x) fails!
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
        ('t_label', REFERENCE),
        ('f_label', REFERENCE),
        ('content', REFERENCE),
        ('t_init', REFERENCE),
    ]

    def as_class(self, cls):
        if not issubclass(cls, Spectrogram):
            raise ValueError

        dct = {}
        var = vars(self)
        for prop, _ in cls.COPY_PROPERTIES:
            if not prop in var:
                raise ValueError
            dct[prop] = var[prop]
        return cls(self, **dct)

    def get_params(self):
        return dict(
            (name, getattr(self, name)) for name, _ in self.COPY_PROPERTIES
        )

    @classmethod
    def _new_with_params(cls, data, params, **kwargs):
        params.update(kwargs)
        """ Implementation detail. """
        obj = cls.__new__(cls, data)
        for key, value in params.iteritems():
            setattr(obj, key, value)
        return obj
    
    def slice(self, y_range, x_range):
        """ Return new spectrogram reduced to the values passed
        as slices. """
        data = super(Spectrogram, self).__getitem__([y_range, x_range])
        params = self.get_params()

        soffset = 0 if x_range.start is None else x_range.start
        eoffset = self.shape[1] if x_range.stop is None else x_range.stop # pylint: disable=E1101
        eoffset -= 1

        fsoffset = 0 if y_range.start is None else y_range.start
        feoffset = self.shape[0] if y_range.stop is None else y_range.stop # pylint: disable=E1101
        
        params.update({
            'time_axis': self.time_axis[
                x_range.start:x_range.stop:x_range.step
            ] - self.time_axis[soffset],
            'freq_axis': self.freq_axis[
                y_range.start:y_range.stop:y_range.step],
            'start': self.start + datetime.timedelta(
                seconds=self.time_axis[soffset]),
            'end': self.start + datetime.timedelta(
                seconds=self.time_axis[eoffset]),
            't_init': self.t_init + self.time_axis[soffset],
        })
        return self._new_with_params(data, params)

    # This accepting arbitrary arguments makes it easier to subclass this.
    def __new__(cls, data, *args, **kwargs):
        return np.asarray(data).view(cls)

    def __init__(self, data, time_axis, freq_axis, start, end, t_init=None,
        t_label="Time", f_label="Frequency", content=""):
        # Because of how object creation works, there is no avoiding
        # unused arguments in this case.
        if t_init is None:
            diff = start - get_day(start)
            t_init = diff.seconds
        self.start = start
        self.end = end

        self.t_label = t_label
        self.f_label = f_label

        self.t_init = t_init

        self.time_axis = time_axis
        self.freq_axis = freq_axis

        self.content = content

    def time_formatter(self, x, pos):
        """ This returns the label for the tick of value x at
        a specified pos on the time axis. """
        # Callback, cannot avoid unused arguments.
        # pylint: disable=W0613
        try:
            return self.format_time(
                self.start + datetime.timedelta(
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

    def plot(self, overlays=[], colorbar=True, min_=None, max_=None, 
        **matplotlib_args):
        # [] as default argument is okay here because it is only read.
        # pylint: disable=W0102,R0914

        data = np.array(self.clip(min_, max_))

        figure = plt.figure(frameon=True)
        axes = figure.add_subplot(111)
        
        params = {
            'origin': 'lower',
            'aspect': 'auto',
        }
        params.update(matplotlib_args)
        im = axes.imshow(np.array(self), **params)
        
        xa = axes.get_xaxis()
        ya = axes.get_yaxis()

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

    def auto_const_bg(self):
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
        return bg.reshape(self.shape[0], 1)

    def subtract_bg(self):
        return self - self.auto_const_bg()

    def randomized_auto_const_bg(self, amount):
        """ Perform constant background subtraction. """
        cols = [randint(0, self.shape[1] - 1) for _ in xrange(amount)]

        # pylint: disable=E1101,E1103
        data = self.astype(to_signed(self.dtype))
        # Subtract average value from every frequency channel.
        tmp = (data - np.average(self, 1).reshape(self.shape[0], 1))
        # Get standard deviation at every point of time.
        # Need to convert because otherwise this class's __getitem__
        # is used which assumes two-dimensionality.
        tmp = tmp[:, cols]
        sdevs = np.asarray(np.std(tmp, 0))

        # Get indices of values with lowest standard deviation.
        cand = sorted(xrange(amount), key=lambda y: sdevs[y])
        # Only consider the best 5 %.
        realcand = cand[:max(1, int(0.05 * len(cand)))]

        # Average the best 5 %
        bg = np.average(self[:, [cols[r] for r in realcand]], 1)

        return bg.reshape(self.shape[0], 1)

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

    def rescale(self, min_=0, max_=1, dtype_=np.dtype('float32')):
        """ Normalize intensities to [min_, max_]. """
        if max_ == min_:
            raise ValueError("Maximum and minimum must be different.")
        if self.max() == self.min():
            raise ValueError("Spectrogram needs to contain distinct values.")
        data = self.astype(dtype_) # pylint: disable=E1101
        return (
            min_ + (max_ - min_) * (data - self.min()) / # pylint: disable=E1101
            (self.max() - self.min()) # pylint: disable=E1101
        )

    def interpolate(self, frequency):
        """ Linearly interpolate intensity at unknown frequency. """
        lfreq, lvalue = None, None
        for freq, value in izip(self.freq_axis, self[:, :]):
            if freq < frequency:
                break
            lfreq, lvalue = freq, value
        else:
            raise ValueError("Frequency not in interpolation range")
        if lfreq is None:
            raise ValueError("Frequency not in interpolation range")
        diff = frequency - freq
        ldiff = lfreq - frequency
        return (ldiff * value + diff * lvalue) / (diff + ldiff)

    @staticmethod
    def merge(items, key=lambda x: x):
        state = {}
        for item in map(iter, items):
            try:
                first = item.next()
            except StopIteration:
                continue
            else:
                state[item] = (first, key(first))

        while state:
            for item, (value, tk) in state.iteritems():
                # Value is biggest.
                if all(tk >= k for it, (v, k)
                    in state.iteritems() if it is not item):
                    yield value
                    break
            try:
                n = item.next()
                state[item] = (n, key(n))
            except StopIteration:
                del state[item]

    def linearize_freqs(self, delta_freq=None):
        # 
        if delta_freq is None:
            delta_freq = (self.freq_axis[:-1] - self.freq_axis[1:])
            # Multiple values at the same frequency are just thrown away
            # in the process of linearizaion
            delta_freq = delta_freq[delta_freq != 0].min()
        nsize = (self.freq_axis.max() - self.freq_axis.min()) / delta_freq + 1
        new = np.zeros((nsize, self.shape[1]), dtype=self.dtype)

        freqs = self.freq_axis - self.freq_axis.min()
        freqs = freqs / delta_freq

        midpoints = np.round((freqs[:-1] + freqs[1:]) / 2)
        fillto = np.concatenate(
            [midpoints, np.round([freqs[-1]])]
        )
        fillfrom = np.concatenate(
            [np.round([freqs[0] + 1]), midpoints]
        )

        for row, from_, to_ in izip(self, fillfrom, fillto):
            new[to_:from_] = row

        vrs = self.get_params()
        vrs.update({
            'freq_axis': np.linspace(
                self.freq_axis.max(), self.freq_axis.min(), nsize
            )
        })

        return self.__class__(new, **vrs)

    def freq_overlap(self, other):
        lower = max(self.freq_axis[-1], other.freq_axis[-1])
        upper = min(self.freq_axis[0], other.freq_axis[0])
        if lower > upper:
            raise ValueError("No overlap.")
        return lower, upper


class LinearTimeSpectrogram(Spectrogram):
    # pylint: disable=E1002
    COPY_PROPERTIES = Spectrogram.COPY_PROPERTIES + [
        ('t_delt', REFERENCE),
    ]
     
    def __init__(self, data, time_axis, freq_axis, start, end,
        t_init, t_delt, t_label="Time", f_label="Frequency",
        content=""):
        super(LinearTimeSpectrogram, self).__init__(
            data, time_axis, freq_axis, start, end, t_init, t_label, f_label,
            content
        )
        self.t_delt = t_delt

    @staticmethod
    def make_array(shape):
        return np.zeros(shape)

    @staticmethod
    def memmap(filename):
        return (
            lambda shape: np.memmap(filename, mode="write", shape=shape)
        )
    
    def resample_time(self, new_delt):
        """ Rescale image so that the difference in time between pixels is
        new_delt seconds. """
        if self.t_delt == new_delt:
            return self
        factor = self.t_delt / float(new_delt)

        # The last data-point does not change!
        new_size = floor((self.shape[1] - 1) * factor + 1) # pylint: disable=E1101
        data = ndimage.zoom(self, (1, new_size / self.shape[1])) # pylint: disable=E1101

        params = self.get_params()
        params.update({
            'time_axis': np.linspace(
                self.time_axis[0],
                self.time_axis[(new_size - 1) * new_delt / self.t_delt],
                new_size
            ),
            't_delt': new_delt,
        })
        return self._new_with_params(data, params)

    @classmethod
    def join_many(cls, spectrograms, mk_arr=None, nonlinear=False,
        maxgap=0, fill=0):
        # XXX: Only load header and load contents of files
        # on demand.
        mask = None

        if mk_arr is None:
            mk_arr = cls.make_array

        specs = sorted(spectrograms, key=lambda x: x.start)

        freqs = specs[0].freq_axis
        if not all(np.array_equal(freqs, sp.freq_axis) for sp in specs):
            raise ValueError("Frequeny channels do not match.")

        # Smallest time-delta becomes the common time-delta.
        min_delt = min(sp.t_delt for sp in specs)

        specs = [sp.resample_time(min_delt) for sp in specs]
        size = sum(sp.shape[1] for sp in specs)

        data = specs[0]
        init = data.t_init
        start_day = data.start

        xs = []
        last = data
        for elem in specs[1:]:
            e_init = (
                SECONDS_PER_DAY * (
                    get_day(elem.start) - get_day(start_day)
                ).days + elem.t_init
            )
            x = int((e_init - last.t_init) / min_delt)
            xs.append(x)
            diff = last.shape[1] - x

            if maxgap is not None and -diff > maxgap:
                raise ValueError("Too large gap.")

            # If we leave out undefined values, we do not want to
            # add values here if x > t_res.
            if nonlinear:
                size -= max(0, diff)
            else:
                size -= diff

            last = elem

        # The non existing element after the last one starts after
        # the last one. Needed to keep implementation below sane.
        xs.append(specs[-1].shape[1])

        # We do that here so the user can pass a memory mapped
        # array if they'd like to.
        arr = mk_arr((data.shape[0], size))
        time_axis = np.zeros((size,))
        sx = 0
        # Amount of pixels left out due to nonlinearity. Needs to be
        # considered for correct time axes.
        sd = 0
        for x, elem in izip(xs, specs):
            diff = x - elem.shape[1]
            if x > elem.shape[1]:
                if nonlinear:
                    x = elem.shape[1]
                else:
                    # If we want to stay linear, fill up the missing
                    # pixels with placeholder zeros.
                    filler = np.zeros((data.shape[0], diff))
                    filler[:] = fill
                    minimum = elem.time_axis[-1]
                    elem = cls._new_with_params(
                        np.concatenate([elem, filler], 1),
                        elem.get_params(),
                        time_axis=np.concatenate(
                            [
                                elem.time_axis,
                                np.linspace(
                                    minimum + min_delt,
                                    minimum + diff * min_delt,
                                    diff
                                )
                            ]
                        )
                    )
            
            arr[:, sx:sx + x] = elem[:, :x]
            if diff > 0:
                if mask is None:
                    mask = np.zeros((data.shape[0], size), dtype=np.uint8)
                mask[:, sx + x - diff:sx + x] = 1
            time_axis[sx:sx + x] = elem.time_axis[:x] + data.t_delt * (sx + sd)
            if nonlinear:
                sd += max(0, diff)
            sx += x
        params = {
            'time_axis': time_axis,
            'freq_axis': data.freq_axis,
            'start': data.start,
            'end': specs[-1].end,
            't_delt': data.t_delt,
            't_init': data.t_init,
            't_label': data.t_label,
            'f_label': data.f_label,
            'content': data.content,
        }
        if mask is not None:
            arr = ma.array(arr, mask=mask)
        if nonlinear:
            del params['t_delt']
            return Spectrogram(arr, **params)
        return LinearTimeSpectrogram(arr, **params)

    def time_to_x(self, time):
        """ Return x-coordinate in spectrogram that corresponds to the
        passed datetime value. """
        # This is impossible for frequencies because that mapping
        # is not injective.
        diff = time - self.start
        diff_s = SECONDS_PER_DAY * diff.days + diff.seconds
        result = diff_s // self.t_delt
        if 0 <= result < self.shape[1]: # pylint: disable=E1101
            return result
        raise ValueError("Out of range.")

    @staticmethod
    def intersect_time(specs):
        delt = min(sp.t_delt for sp in specs)
        start = max(sp.t_init for sp in specs)

        # XXX: Could do without resampling by using
        # sp.t_init below, not sure if good idea.
        specs = [sp.resample_time(delt) for sp in specs]
        cut = [sp[:, (start - sp.t_init) / delt:] for sp in specs]

        length = min(sp.shape[1] for sp in cut)
        return [sp[:, :length] for sp in cut]

    @classmethod
    def combine_frequencies(cls, specs):
        specs = cls.intersect_time(specs)

        one = specs[0]

        dtype_ = max(sp.dtype for sp in specs)
        fsize = sum(sp.shape[0] for sp in specs)

        new = np.zeros((fsize, one.shape[1]), dtype=dtype_)

        freq_axis = np.zeros((fsize,))


        for n, (data, row) in enumerate(cls.merge(
            [
                [(sp, n) for n in xrange(sp.shape[0])] for sp in specs
            ],
            key=lambda x: x[0].freq_axis[x[1]]
        )):
            new[n, :] = data[row, :]
            freq_axis[n] = data.freq_axis[row]
        params = {
            'time_axis': one.time_axis, # Should be equal
            'freq_axis': freq_axis,
            'start': one.start,
            'end': one.end,
            't_delt': one.t_delt,
            't_init': one.t_init,
            't_label': one.t_label,
            'f_label': one.f_label,
            'content': one.content,
        }
        # XXX
        return cls._new_with_params(new, params)
