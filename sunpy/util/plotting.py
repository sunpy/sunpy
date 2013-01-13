# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 15:34:49 2013

@author: Stuart Mumford

Some Independant plotting tools, mainly animation UI based.
"""
__author__ = "Stuart Mumford"
__email__ = "stuartmumford@physics.org"

import matplotlib.animation as animation
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mpl_toolkits.axes_grid1.axes_size as Size

class ControlFuncAnimation(animation.FuncAnimation):
    """ This is a slight modification to the animation class to allow pausing
    starting and stopping."""
    def __init__(self, fig, func, frames=None, init_func=None, fargs=None,
            save_count=None, auto_start=True, **kwargs):
        self.fig = fig #This should be done.
        animation.FuncAnimation.__init__(self, fig, func, frames=frames,
                                         init_func=init_func, fargs=fargs,
                                         save_count=save_count, **kwargs)
        if not auto_start:
            self._fig.canvas.mpl_disconnect(self._first_draw_id)
            self._first_draw_id = None
    
    def _start(self, *args):
        if self.event_source is None:
            self.event_source = self.fig.canvas.new_timer()
            self.event_source.interval = self._interval
        animation.FuncAnimation._start(self)
    
    def _stop(self, *args):
        if self.event_source:
            animation.FuncAnimation._stop(self, *args)

def add_controls(axes=None):
    """ Adds Start/Stop controls to an axes having been given a animation 
    instance. """
    
    #If No axes specified use current axes.
    if not axes:
        axes = plt.gca()
    fig = axes.get_figure()
    
    #Split up the current axes so there is space for a start and a stop button
    divider = make_axes_locatable(axes)
    pad = 0.1 # Padding between axes
    pad_size = Size.Fraction(pad, Size.AxesX(axes))

    #Define size of usefult axes cells, 50% each in x 20% for buttons in y.
    xsize = Size.Fraction((1.-2.*pad)/2., Size.AxesX(axes))
    ysize = Size.Fraction((1.-2.*pad)/5., Size.AxesY(axes))

    #Set up grid, 3x3 with cells for padding.
    divider.set_horizontal([xsize, pad_size, xsize])
    divider.set_vertical([ysize, pad_size, Size.AxesY(axes)])

    #Main figure spans all horiz and is in the top (2) in vert.
    axes.set_axes_locator(divider.new_locator(0, 2, nx1=-1))
    
    #Add two axes for buttons and make them 50/50 spilt at the bottom.
    bax1 = fig.add_axes((0.,0.,1.,1.))
    locator = divider.new_locator(nx=0, ny=0)
    bax1.set_axes_locator(locator)
    bax2 = fig.add_axes((0.,0.,0.8,1.))
    locator = divider.new_locator(nx=2, ny=0)
    bax2.set_axes_locator(locator)
    
    return axes, bax1, bax2