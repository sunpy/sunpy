# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101, E0611
'''
Main class for representing cubes - 3D sets of continuous data by time and/or
wavelength
'''
# NOTE: This module uses version 1.02 of "Time coordinates in FITS" by
# Rots et al, available at http://hea-www.cfa.harvard.edu/~arots/TimeWCS/
# This draft standard may change.

# standard libraries
import datetime

# external libraries
import numpy as np
import matplotlib.pyplot as plt
import astropy.nddata
import astropy.units as u
from astropy.units import sday  # sidereal day

# Sunpy modules
from sunpy.map import GenericMap, MapMeta
from sunpy.visualization.imageanimator import ImageAnimator
from sunpy.lightcurve import LightCurve
from sunpy.spectra.spectrum import Spectrum
from sunpy.spectra.spectrogram import Spectrogram
from sunpy.cube import cube_utils as cu

__all__ = ['Cube']


class Cube(astropy.nddata.NDData):
    ''' Class representing spectral cubes.

        Attributes
        ----------
        data: numpy ndarray
            The spectral cube holding the actual data in this object. The axes'
            priorities are time, spectral, celestial. This means that if
            present, each of these axis will take precedence over the others.
            For example, in an x, y, t cube the order would be (t,x,y) and in a
            lambda, t, y cube the order will be (t, lambda, y).

        axes_wcs: sunpy.wcs.wcs.WCS object
            The WCS object containing the axes' information

        meta: dict
            Header containing the wavelength-specific metadata as well as the
            whole-file metadata
    '''

    def __init__(self, data, wcs, meta=None, **kwargs):
        data, wcs = cu.orient(data, wcs)
        astropy.nddata.NDData.__init__(self, data=data,
                                       meta=meta,
                                       **kwargs)
        self.axes_wcs = wcs
        # We don't send this to NDData because it's not
        # supported as of astropy 0.3.2. Eventually we will.
        # Also it's called axes_wcs because wcs belongs to astropy.nddata and
        # that messes up slicing.

    def plot_wavelength_slice(self, offset, axes=None,
                              style='imshow', **kwargs):
        '''
        Plots an x-y graph at a certain specified wavelength onto the current
        axes. Keyword arguments are passed on to matplotlib.

        Parameters
        ----------
        offset: int or float
            The offset from the primary wavelength to plot. If it's an int it
            will plot the nth wavelength from the primary; if it's a float then
            it will plot the closest wavelength. If the offset is out of range,
            it will plot the primary wavelength (offset 0)

        axes: matplotlib.axes or None:
            The axes to plot onto. If None the current axes will be used.

        style: 'imshow' or 'pcolormesh'
            The style of plot to be used. Default is 'imshow'
        '''
        if axes is None:
            axes = plt.gca()

        data = self._choose_wavelength_slice(offset)
        if data is None:
            data = self._choose_wavelength_slice(0)

        if style is 'imshow':
            plot = axes.imshow(data, **kwargs)
        elif style is 'pcolormesh':
            plot = axes.pcolormesh(data, **kwargs)

        return plot

    def plot_x_slice(self, offset, axes=None,
                     style='imshow', **kwargs):
        '''
        Plots an x-y graph at a certain specified wavelength onto the current
        axes. Keyword arguments are passed on to matplotlib.

        Parameters
        ----------
        offset: int or float
            The offset from the initial x value to plot. If it's an int it
            will plot slice n from the start; if it's a float then
            it will plot the closest x-distance. If the offset is out of range,
            it will plot the primary wavelength (offset 0)

        axes: matplotlib.axes or None:
            The axes to plot onto. If None the current axes will be used.

        style: 'imshow' or 'pcolormesh'
            The style of plot to be used. Default is 'imshow'
        '''
        if axes is None:
            axes = plt.gca()

        data = self._choose_x_slice(offset)
        if data is None:
            data = self._choose_x_slice(0)

        if style is 'imshow':
            plot = axes.imshow(data, **kwargs)
        elif style is 'pcolormesh':
            plot = axes.pcolormesh(data, **kwargs)

        return plot

    def animate(self, *args, **kwargs):
        '''Plots an interactive visualization of this cube with a slider
        controlling the wavelength axis.
        Parameters other than data are passed to ImageAnimator, which in turn
        passes them to imshow.'''
        i = ImageAnimator(data=self.data, *args, **kwargs)
        return i

    def _choose_wavelength_slice(self, offset):
        '''Retrieves an x-y slice at a wavelength specified by the cube's
        primary wavelength plus the given offset.

        Parameters
        ----------
        offset: int or astropy quantity
            Offset from the cube's primary wavelength. If the value is an int,
            then it returns that slice. Otherwise, it will return the nearest
            wavelength to the one specified.
        '''
        if 'WAVE' not in self.axes_wcs.wcs.ctype:
            raise cu.CubeError(2, "Spectral dimension not present")

        axis = 1 if self.axes_wcs.wcs.ctype[-1] in ['TIME', 'UTC'] else 0
        arr = None
        length = self.data.shape[axis]
        if isinstance(offset, int) and offset >= 0 and offset < length:
            arr = self.data.take(offset, axis=axis)

        if isinstance(offset, u.Quantity):
            delta = self.axes_wcs.wcs.cdelt[-1 - axis] * u.m
            wloffset = offset.to(u.m) / delta
            wloffset = int(wloffset)
            if wloffset >= 0 and wloffset < self.data.shape[axis]:
                arr = self.data.take(wloffset, axis=axis)

        return arr

    def _choose_x_slice(self, offset):
        '''
        Retrieves a lambda-y slice at an x coordinate specified by the cube's
        primary wavelength plus the given offset.

        Parameters
        ----------
        offset: int or astropy quantity
            Offset from the cube's initial x. If the value is an int,
            then it returns that slice. Otherwise, it will return the nearest
            wavelength to the one specified.
        '''
        arr = None
        axis = 1 if self.axes_wcs.wcs.ctype[-2] != 'WAVE' else 2
        length = self.data.shape[axis]
        if isinstance(offset, int) and offset >= 0 and offset < length:
            arr = self.data.take(offset, axis=axis)

        if isinstance(offset, u.Quantity):
            unit = self.axes_wcs.wcs.cunit[-1 - axis]
            delta = self.axes_wcs.wcs.cdelt[-1 - axis] * unit
            wloffset = offset.to(unit) / delta
            wloffset = int(wloffset)
            if wloffset >= 0 and wloffset < self.data.shape[axis]:
                arr = self.data.take(wloffset, axis=axis)

        return arr

    def slice_to_map(self, chunk, *args, **kwargs):
        # TODO: implement slice-by-float functionality
        '''
        Converts a given frequency chunk to a SunPy Map. Extra parameters are
        passed on to Map.

        Parameters
        ----------
        chunk: int or float or (int, int) or (float, float)
            The piece of the cube to convert to a map. If it's a single number,
            then it will return that single-slice map, otherwise it will
            aggregate the given range.
        '''
        if self.axes_wcs.wcs.ctype[-2] == 'WAVE':
            error = "Cannot construct a map with only one spatial dimension"
            raise cu.CubeError(3, error)

        if isinstance(chunk, tuple):
            maparray = self.data[chunk[0]:chunk[1], :, :].sum(0)
        else:
            maparray = self.data[chunk, :, :]
        mapheader = MapMeta(self.meta)
        gmap = GenericMap(data=maparray, header=mapheader, *args, **kwargs)
        return gmap

    def slice_to_lightcurve(self, wavelength, y_coord=None):
        '''
        For a time-lambda-y cube, returns a lightcurve with curves at the
        specified wavelength and given y-coordinate. If no y is given, all of
        them will be used (meaning the lightcurve object could contain more
        than one timecurve.)
        Parameters
        ----------
        wavelength: int or float
            The wavelength to take the y-coordinates from
        y_coord: int, optional
            The y-coordinate to take the lightcurve from.
        '''
        if self.axes_wcs.wcs.ctype[-1] not in ['TIME', 'UTC']:
            raise cu.CubeError(1,
                               'Cannot create a lightcurve with no time axis')
        if self.axes_wcs.wcs.ctype[-2] != 'WAVE':
            raise cu.CubeError(2, 'A spectral axis is needed in a lightcurve')

        data = self._choose_wavelength_slice(wavelength)
        if y_coord is not None:
            data = data[:, y_coord]
        return LightCurve(data=data, meta=self.meta)

    def slice_to_spectrum(self, fst_coord, snd_coord):
        '''
        For a cube containing a spectral dimension, returns a sunpy spectrum.
        The given coordinates represent which values to take. If they are None,
        then the corresponding axis is summed.
        Parameters
        ----------
        fst_coord: int or None
            The first coordinate to pick. Keep in mind that depending on the
            cube, this may be in the first or second axis. If None, the whole
            axis will be taken and its values summed.
        snd_coord: int or None
            The second coordinate to pick. This will always correspond to the
            third axis. If None, the whole axis will be taken and its values
            summed.
        '''
        if 'WAVE' not in self.axes_wcs.wcs.ctype:
            raise cu.CubeError(2, 'Spectral axis needed to create a spectrum')
        axis = 0 if self.axes_wcs.wcs.ctype[-1] == 'WAVE' else 1

        if axis == 0:
            data = self.data[:, fst_coord, snd_coord]
        else:
            data = self.data[fst_coord, :, snd_coord]

        if snd_coord is None:
            data = data.sum(axis=2)
        if fst_coord is None:
            sumaxis = 1 if axis == 0 else 0
            data = data.sum(axis=sumaxis)

        freq_axis = self.freq_axis()
        return Spectrum(np.array(data), np.array(freq_axis))

    def slice_to_spectrogram(self, y_coord, **kwargs):
        # TODO: make this not take only an int.
        '''
        For a time-lambda-y cube, given a y-coordinate, returns a sunpy
        spectrogram. Keyword arguments are passed on to Spectrogram's __init__.
        Parameters
        ----------
        y_coord: int
            The y-coordinate to pick when converting to a spectrogram.
        '''
        if self.axes_wcs.wcs.ctype[-1] not in ['TIME', 'UTC']:
            raise cu.CubeError(1,
                               'Cannot create a spectrogram with no time axis')
        if self.axes_wcs.wcs.ctype[-2] != 'WAVE':
            raise cu.CubeError(2, 'A spectral axis is needed in a spectrogram')
        data = self.data[:, :, y_coord]
        time_axis = self.time_axis()
        freq_axis = self.freq_axis()

        if 'DATE_OBS'in self.meta:
            tformat = '%Y-%m-%dT%H:%M:%S.%f'
            start = datetime.datetime.strptime(self.meta['DATE_OBS'], tformat)
        else:
            start = datetime.datetime(1, 1, 1)

        if 'DATE_END' in self.meta:
            tformat = '%Y-%m-%dT%H:%M:%S.%f'
            end = datetime.datetime.strptime(self.meta['DATE_END'], tformat)
        else:
            dif = time_axis[-1] - time_axis[0]
            unit = self.axes_wcs.wcs.cunit[-1]
            dif = dif * u.Unit(unit)
            days = dif.to(sday)
            lapse = datetime.timedelta(days.value)
            end = start + lapse
        return Spectrogram(data=data, time_axis=time_axis, freq_axis=freq_axis,
                           start=start, end=end, **kwargs)

    def time_axis(self):
        '''
        Returns a numpy array containing the time values for the cube's time
        dimension.
        '''
        if self.axes_wcs.wcs.ctype[-1] not in ['TIME', 'UTC']:
            raise cu.CubeError(1, 'No time axis present')
        delta = self.axes_wcs.wcs.cdelt[-1]
        crpix = self.axes_wcs.wcs.crpix[-1]
        crval = self.axes_wcs.wcs.crval[-1]
        start = crval - crpix * delta
        stop = start + len(self.data) * delta
        return np.arange(start, stop, delta)

    def freq_axis(self):
        '''
        Returns a numpy array containing the frequency values for the cube's
        spectral dimension.
        '''
        if 'WAVE' not in self.axes_wcs.wcs.ctype:
            raise cu.CubeError(2,
                               'No energy (wavelength, frequency) axis found')
        axis = 0 if self.axes_wcs.wcs.ctype[-1] == 'WAVE' else 1
        delta = self.axes_wcs.wcs.cdelt[-1 - axis]
        crpix = self.axes_wcs.wcs.crpix[-1 - axis]
        crval = self.axes_wcs.wcs.crval[-1 - axis]
        start = crval - crpix * delta
        stop = start + self.data.shape[axis] * delta
        return np.arange(start, stop, delta)

    def __getitem__(self, item):
        if item is None or (isinstance(item, tuple) and None in item):
            raise IndexError("None indices not supported")
        axes = self.axes_wcs.wcs.ctype
        slice_to_map = (axes[-2] != 'WAVE' and
                        (isinstance(item, int) or
                         cu.iter_isinstance(item, int, slice) or
                         cu.iter_isinstance(item, int, slice, slice)))
        slice_to_spectrum = (((isinstance(item, int) or
                               cu.iter_isinstance(item, int, slice) or
                               cu.iter_isinstance(item, int, slice, slice) or
                               cu.iter_isinstance(item, int, slice, int))
                              and axes[-2] == 'WAVE')
                             or (axes[-1] == 'WAVE' and
                                 cu.iter_isinstance(item, slice, int, int)))
        slice_to_spectrogram = (cu.iter_isinstance(item, slice, slice, int)
                                and axes[-2] == 'WAVE')
        slice_to_lightcurve = (axes[-2] == 'WAVE' and
                               (cu.iter_isinstance(item, slice, int, int)
                                or cu.iter_isinstance(item, slice, int)
                                or cu.iter_isinstance(item,
                                                      slice, int, slice)))
        stay_as_cube = (isinstance(item, slice) or
                        (isinstance(item, tuple) and
                         not any(isinstance(i, int) for i in item)))

        reducedcube = cu.reduce_dim(self, 0, slice(None, None, None))
        if isinstance(item, tuple):
            for i in range(len(item)):
                if isinstance(item[i], slice):
                    reducedcube = cu.reduce_dim(reducedcube, i, item[i])

        if slice_to_map:
            return cu.handle_slice_to_map(reducedcube, item)
        elif slice_to_spectrum:
            return cu.handle_slice_to_spectrum(reducedcube, item)
        elif slice_to_spectrogram:
            return reducedcube.slice_to_spectrogram(item[2])
        elif slice_to_lightcurve:
            return cu.handle_slice_to_lightcurve(reducedcube, item)
        elif stay_as_cube:
            return reducedcube
        else:
            return self.data[item]
