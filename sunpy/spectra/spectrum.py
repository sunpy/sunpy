# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

import numpy as np
from sunpy.time import parse_time

from matplotlib import pyplot as plt

class Spectrum(np.ndarray):
    """ Class representing a general spectrum.
    
    Parameters
    ----------
    flux : np.ndarray
        one-dimensional array with the flux intensity in a particular channel.
    channel : 2d np.ndarray
        two-dimensional array with the channel boundaries [lo, hi]
    width : np.ndarray
        bin width of the channels (allowed to be zero).
    error : np.ndarray
        error in flux in each channel
    units : 2d str
        units for the channel (e.g. energy, wavelength) and flux (e.g. photons/cm^2/s^1)
    name : str
        name of the spectrum
    """

    def __new__(cls, data, *args, **kwargs):
        return np.asarray(data).view(cls)

    def __init__(self, channel, flux, width=None, error=None, name='spectrum', units=None, date=None):      
        self.flux = flux
        self.name = name
        self.units = units
        self.width = width
        if date is not None: self.date = parse_time(date) 
        else: self.date = None

        if units == None:
            self.units = ['frequency [MHz]', 'solar flux units [sfu]']
        
        if np.array(channel.shape).shape[0] == 1:
            
            if channel.shape != flux.shape: 
                raise ValueError("channel and flux don't match shapes!")

            if width != None:
                self.width = width
                channel1 = channel
                channel2 = channel + width
                self.channel = np.array([channel, channel + width])

            if width == None:
                self.channel = np.array([channel, channel])
        
        if np.array(channel.shape).shape[0] == 2:
            if channel.shape[1] != flux.shape[0]:
                raise ValueError("channel and flux don't match shapes!")
            self.channel = channel
            self.width = channel[1,:] - channel[0,:]

    def plot(self, figure=None, overlays=[], log=None, **matplotlib_args):
        """
        Plot spectrum onto figure.
        
        Parameters
        ----------
        figure : matplotlib.figure.Figure
            Figure to plot the spectrogram on. If None, new Figure is created.
        overlays : list
            List of overlays (functions that receive figure and axes and return
            new ones) to be applied after drawing.
        """
        # [] as default argument is okay here because it is only read.
        # pylint: disable=W0102,R0914
        if figure is None:
            figure = plt.figure(frameon=True)
        axes = figure.add_subplot(111)

        params = {}
        params.update(matplotlib_args)

        # todo: add code to detect whether log scaling is necessary
        axes.set_yscale('log')
        axes.set_xscale('log')
   
        if self.width is None:
            axes.plot(self.channel[0,:], self.flux, **params)
        if self.width is not None:
            axes.step(self.channel[0,:], self.flux, where='post', **params)

        if ((self.name != None) & (self.date != None)):
            axes.set_title("%s %s" % (self.name, self.date))
        if self.name != None:
            axes.set_title("%s" % self.name)
        if self.date != None:
            axes.set_title("%s" % self.date)

        axes.set_xlabel(self.units[0])
        axes.set_ylabel(self.units[1])
            
        #for overlay in overlays:
        #    figure, axes = overlay(figure, axes)
        return figure
        
    def show(self, figure=None, overlays=None, **matplotlib_args):
        """Displays map on screen. Arguments are same as plot()."""
        self.plot(figure, overlays).show()
