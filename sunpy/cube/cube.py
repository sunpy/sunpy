# -*- coding: utf-8 -*-
# Author: Mateo Inchaurrandieta <mateo.inchaurrandieta@gmail.com>
# pylint: disable=E1101, E0611
"""
Main class for representing cubes - 3D sets of continuous data by time and/or
wavelength
"""
# NOTE: This module uses version 1.02 of "Time coordinates in FITS" by
# Rots et al, available at http://hea-www.cfa.harvard.edu/~arots/TimeWCS/
# This draft standard may change.

# standard libraries
import datetime

# external libraries
import numpy as np
import matplotlib.pyplot as plt
import astropy.nddata
from astropy import units as u
from astropy.units import sday  # sidereal day

# Sunpy modules
from sunpy.map import GenericMap, MapMeta
from sunpy.visualization.imageanimator import ImageAnimator
from sunpy.lightcurve import LightCurve
from sunpy.spectra.spectrum import Spectrum
from sunpy.spectra.spectrogram import Spectrogram
from sunpy.spectra.spectral_cube import SpectralCube
from sunpy.cube import cube_utils as cu
from sunpy.wcs import wcs_util as wu

__all__ = ['Cube']


class Cube(astropy.nddata.NDDataArray):
    """
    Class representing cubes.
    Extra arguments are passed on to NDDataArray's init.

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
    """

    def __init__(self, data, wcs, **kwargs):
        data, wcs = cu.orient(data, wcs)
        astropy.nddata.NDDataArray.__init__(self, data=data, **kwargs)
        self.axes_wcs = wcs
        # We don't send this to NDDataArray because it's not
        # supported as of astropy 1.0. Eventually we will.
        # Also it's called axes_wcs because wcs belongs to astropy.nddata and
        # that messes up slicing.

    def plot_wavelength_slice(self, offset, axes=None,
                              style='imshow', **kwargs):
        """
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
        """
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
        """
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
        """
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
        """
        Plots an interactive visualization of this cube with a slider
        controlling the wavelength axis.
        Parameters other than data are passed to ImageAnimator, which in turn
        passes them to imshow.
        """
        i = ImageAnimator(data=self.data, *args, **kwargs)
        return i

    def _choose_wavelength_slice(self, offset):
        """
        Retrieves an x-y slice at a wavelength specified by the cube's
        primary wavelength plus the given offset.

        Parameters
        ----------
        offset: int or astropy quantity
            Offset from the cube's primary wavelength. If the value is an int,
            then it returns that slice. Otherwise, it will return the nearest
            wavelength to the one specified.
        """
        if 'WAVE' not in self.axes_wcs.wcs.ctype:
            raise cu.CubeError(2, "Spectral dimension not present")
        if self.data.ndim == 4:
            raise cu.CubeError(4, "Can only work with 3D cubes")

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
        """
        Retrieves a lambda-y slice at an x coordinate specified by the cube's
        primary wavelength plus the given offset.

        Parameters
        ----------
        offset: int or astropy quantity
            Offset from the cube's initial x. If the value is an int,
            then it returns that slice. Otherwise, it will return the nearest
            wavelength to the one specified.
        """
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

    def slice_to_map(self, chunk, snd_dim=None, *args, **kwargs):
        """
        Converts a given frequency chunk to a SunPy Map. Extra parameters are
        passed on to Map.

        Parameters
        ----------
        chunk: int or astropy quantity or tuple
            The piece of the cube to convert to a map. If it's a single number,
            then it will return that single-slice map, otherwise it will
            aggregate the given range. Depending on the cube, this may
            correspond to a time or an energy dimension
        snd_dim: int or astropy quantity or tuple, optional
            Only used for hypercubes, the wavelength to choose from; works in
            the same way as chunk.
        """
        if self.axes_wcs.wcs.ctype[-2] == 'WAVE' and self.data.ndim == 3:
            error = "Cannot construct a map with only one spatial dimension"
            raise cu.CubeError(3, error)

        def pixelize(coord, axis):
            '''shorthand for convert_point'''
            unit = coord.unit if isinstance(coord, u.Quantity) else None
            return cu.convert_point(coord, unit, self.axes_wcs, axis)
        if isinstance(chunk, tuple):
            item = slice(pixelize(chunk[0], 0), pixelize(chunk[1], 0), None)
            maparray = self.data[item].sum(0)
        else:
            maparray = self.data[pixelize(chunk, 0)]

        if self.data.ndim == 4:
            if snd_dim is None:
                error = "snd_dim must be given when slicing hypercubes"
                raise cu.CubeError(4, error)

            if isinstance(snd_dim, tuple):
                item = slice(pixelize(snd_dim[0], 1),
                             pixelize(snd_dim[1], 1), None)
                maparray = maparray[item].sum(0)
            else:
                maparray = maparray[pixelize(snd_dim, 1)]

        mapheader = MapMeta(self.meta)
        gmap = GenericMap(data=maparray, header=mapheader, *args, **kwargs)
        return gmap

    def slice_to_lightcurve(self, wavelength, y_coord=None, x_coord=None):
        """
        For a time-lambda-y cube, returns a lightcurve with curves at the
        specified wavelength and given y-coordinate. If no y is given, all of
        them will be used (meaning the lightcurve object could contain more
        than one timecurve.)
        Parameters
        ----------
        wavelength: int or astropy quantity
            The wavelength to take the y-coordinates from
        y_coord: int or astropy quantity, optional
            The y-coordinate to take the lightcurve from.
        x_coord: int or astropy quantity, optional
            In the case of hypercubes, specify an extra celestial coordinate.
        """
        if self.axes_wcs.wcs.ctype[-1] not in ['TIME', 'UTC']:
            raise cu.CubeError(1,
                               'Cannot create a lightcurve with no time axis')
        if self.axes_wcs.wcs.ctype[-2] != 'WAVE':
            raise cu.CubeError(2, 'A spectral axis is needed in a lightcurve')

        def pixelize(coord, axis):
            '''shorthand for convert_point'''
            unit = coord.unit if isinstance(coord, u.Quantity) else None
            return cu.convert_point(coord, unit, self.axes_wcs, axis)

        if self.data.ndim == 3:
            data = self._choose_wavelength_slice(wavelength)
            if y_coord is not None:
                data = data[:, pixelize(y_coord, 1)]
        else:
            if y_coord is None and x_coord is None:
                raise cu.CubeError(4, "At least one coordinate must be given")
            if y_coord is None:
                y_coord = slice(None, None, None)
            else:
                y_coord = pixelize(y_coord, 2)
            if x_coord is None:
                x_coord = slice(None, None, None)
            else:
                x_coord = pixelize(x_coord, 3)
            item = (slice(None, None, None), wavelength, y_coord, x_coord)
            data = self.data[item]

        return LightCurve(data=data, meta=self.meta)

    def slice_to_spectrum(self, *coords):
        """
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
        """
        if 'WAVE' not in self.axes_wcs.wcs.ctype:
            raise cu.CubeError(2, 'Spectral axis needed to create a spectrum')
        axis = 0 if self.axes_wcs.wcs.ctype[-1] == 'WAVE' else 1

        def pixelize(coord):
            '''shorthand for convert_point'''
            unit = coord.unit if isinstance(coord, u.Quantity) else None
            return cu.convert_point(coord, unit, self.axes_wcs, axis)
        pixels = [pixelize(coord) for coord in coords]
        item = range(len(pixels))
        if axis == 0:
            item[1:] = pixels
            item[0] = slice(None, None, None)
            item = [slice(None, None, None) if i is None else i for i in item]
        else:
            item[0] = pixels[0]
            item[1] = slice(None, None, None)
            item[2:] = pixels[1:]
            item = [slice(None, None, None) if i is None else i for i in item]

        data = self.data[item]
        for i in range(len(pixels)):
            if pixels[i] is None:
                if i == 0:
                    sumaxis = 1 if axis == 0 else 0
                else:
                    sumaxis = 1 if i == 2 else i
                data = data.sum(axis=sumaxis)

        freq_axis, cunit = self.freq_axis()
        return Spectrum(np.array(data), np.array(freq_axis), cunit)

    def slice_to_spectrogram(self, y_coord, x_coord=None, **kwargs):
        """
        For a time-lambda-y cube, given a y-coordinate, returns a sunpy
        spectrogram. Keyword arguments are passed on to Spectrogram's __init__.

        Parameters
        ----------
        y_coord: int
            The y-coordinate to pick when converting to a spectrogram.
        x_coord: int
            The x-coordinate to pick. This is only used for hypercubes.
        """
        if self.axes_wcs.wcs.ctype[-1] not in ['TIME', 'UTC']:
            raise cu.CubeError(1,
                               'Cannot create a spectrogram with no time axis')
        if self.axes_wcs.wcs.ctype[-2] != 'WAVE':
            raise cu.CubeError(2, 'A spectral axis is needed in a spectrogram')

        def pixelize(coord, axis):
            '''shorthand for convert_point'''
            unit = coord.unit if isinstance(coord, u.Quantity) else None
            return cu.convert_point(coord, unit, self.axes_wcs, axis)
        if self.data.ndim == 3:
            data = self.data[:, :, pixelize(y_coord, 2)]
        else:
            if x_coord is None:
                raise cu.CubeError(4, 'An x-coordinate is needed for 4D cubes')
            data = self.data[:, :, pixelize(y_coord, 2), pixelize(x_coord, 3)]
        time_axis = self.time_axis()[0]
        freq_axis = self.freq_axis()[0]

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

    def slice_to_cube(self, axis, chunk, **kwargs):
        """
        For a hypercube, return a 3-D cube that has been cut along the given
        axis and with data corresponding to the given chunk.

        Parameters
        ----------
        axis: int
            The axis to cut from the hypercube
        chunk: int, astropy Quantity or tuple:
            The data to take from the axis
        """
        if self.data.ndim == 3:
            raise cu.CubeError(4, 'Can only slice a hypercube into a cube')

        item = [slice(None, None, None) for _ in range(4)]
        if isinstance(chunk, tuple):
            if cu.iter_isinstance(chunk, (u.Quantity, u.Quantity)):
                pixel0 = cu.convert_point(chunk[0].value, chunk[0].unit,
                                          self.axes_wcs, axis)
                pixel1 = cu.convert_point(chunk[1].value, chunk[1].unit,
                                          self.axes_wcs, axis)
                item[axis] = slice(pixel0, pixel1, None)
            elif cu.iter_isinstance((chunk, int, int)):
                item[axis] = slice(chunk[0], chunk[1], None)
            else:
                raise cu.CubeError(5, "Parameters must be of the same type")
            newdata = self.data[item].sum(axis)
        else:
            unit = chunk.unit if isinstance(chunk, u.Quantity) else None
            pixel = cu.convert_point(chunk, unit, self.axes_wcs, axis)
            item[axis] = pixel
            newdata = self.data[item]
        wcs_indices = [0, 1, 2, 3]
        wcs_indices.remove(3 - axis)
        newwcs = wu.reindex_wcs(self.axes_wcs, np.array(wcs_indices))
        if axis == 2 or axis == 3:
            newwcs = wu.add_celestial_axis(newwcs)
            newwcs.was_augmented = True
        cube = Cube(newdata, newwcs, meta=self.meta, **kwargs)
        return cube

    def convert_to_spectral_cube(self):
        """
        Converts this cube into a SpectralCube. It will only work if the cube
        has exactly three dimensions and one of those is a spectral axis.
        """
        if self.data.ndim == 4:
            raise cu.CubeError(4, "Too many dimensions: Can only convert a " +
                               "3D cube. Slice the cube before converting")
        if 'WAVE' not in self.axes_wcs.wcs.ctype:
            raise cu.CubeError(2, 'Spectral axis needed to create a spectrum')
        axis = 0 if self.axes_wcs.wcs.ctype[-1] == 'WAVE' else 1
        coordaxes = [1, 2] if axis == 0 else [0, 2]  # Non-spectral axes
        newwcs = wu.reindex_wcs(self.axes_wcs, coordaxes)
        time_or_x_size = self.data.shape[coordaxes[0]]
        y_size = self.data.shape[coordaxes[1]]
        spectra = np.empty((time_or_x_size, y_size), dtype=Spectrum)
        for i in range(time_or_x_size):
            for j in range(y_size):
                spectra[i][j] = self.slice_to_spectrum(i, j)
        return SpectralCube(spectra, newwcs, self.meta)

    def time_axis(self):
        """
        Returns a numpy array containing the time values for the cube's time
        dimension, as well as the unit used.
        """
        if self.axes_wcs.wcs.ctype[-1] not in ['TIME', 'UTC']:
            raise cu.CubeError(1, 'No time axis present')
        delta = self.axes_wcs.wcs.cdelt[-1]
        crpix = self.axes_wcs.wcs.crpix[-1]
        crval = self.axes_wcs.wcs.crval[-1]
        start = crval - crpix * delta
        stop = start + len(self.data) * delta
        cunit = u.Unit(self.axes_wcs.wcs.cunit[-1])
        return np.arange(start, stop, delta), cunit

    def freq_axis(self):
        """
        Returns a numpy array containing the frequency values for the cube's
        spectral dimension, as well as the axis's unit.
        """
        if 'WAVE' not in self.axes_wcs.wcs.ctype:
            raise cu.CubeError(2,
                               'No energy (wavelength, frequency) axis found')
        axis = 0 if self.axes_wcs.wcs.ctype[-1] == 'WAVE' else 1
        delta = self.axes_wcs.wcs.cdelt[-1 - axis]
        crpix = self.axes_wcs.wcs.crpix[-1 - axis]
        crval = self.axes_wcs.wcs.crval[-1 - axis]
        start = crval - crpix * delta
        stop = start + self.data.shape[axis] * delta
        cunit = u.Unit(self.axes_wcs.wcs.cunit[-1 - axis])
        return np.arange(start, stop, delta), cunit

    def _array_is_aligned(self):
        """
        Returns whether the wcs system and the array are well-aligned.
        """
        rot_matrix = self.axes_wcs.wcs.pc
        return np.allclose(rot_matrix, np.eye(self.axes_wcs.wcs.naxis))

    def __getitem__(self, item):
        if item is None or (isinstance(item, tuple) and None in item):
            raise IndexError("None indices not supported")
        if not self._array_is_aligned():
            raise cu.CubeError(6, "Slicing on unaligned wcs-array systems " +
                               "not supported at the moment")
        pixels = cu.pixelize_slice(item, self.axes_wcs)
        if self.data.ndim == 3:
            return cu.getitem_3d(self, pixels)
        else:
            return cu.getitem_4d(self, pixels)
