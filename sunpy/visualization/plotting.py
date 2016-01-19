# -*- coding: utf-8 -*-
"""
Some Independent plotting tools, mainly animation UI based.
"""
__author__ = "Stuart Mumford"
__email__ = "stuartmumford@physics.org"
__all__ = ['ControlFuncAnimation', 'add_controls']

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
from mpl_toolkits.axes_grid1 import make_axes_locatable
import mpl_toolkits.axes_grid1.axes_size as Size

class ControlFuncAnimation(animation.FuncAnimation):
    """ This is a slight modification to the animation class to allow pausing
    starting and stopping.

    .. todo::
        improve documentation
    """
    def __init__(self, fig, func, frames=None, init_func=None, fargs=None,
            save_count=None, auto_start=True, **kwargs):
        self.fig = fig #This should be done.
        animation.FuncAnimation.__init__(self, fig, func, frames=frames,
                                         init_func=init_func, fargs=fargs,
                                         save_count=save_count, **kwargs)

        self._started = False #Set to false _start will start animation
        if not auto_start:
            self._fig.canvas.mpl_disconnect(self._first_draw_id)
            self._first_draw_id = None

    def _start(self, *args):
        if not self._started:
            if self.event_source is None:
                self.event_source = self.fig.canvas.new_timer()
                self.event_source.interval = self._interval
            animation.FuncAnimation._start(self)
            self._started = True

    def _stop(self, *args):
        if self.event_source:
            animation.FuncAnimation._stop(self, *args)
        self._started = False

def add_controls(axes=None, slider=False):
    """ Adds Start/Stop controls to an axes having been given a animation
    instance.

    Parameters
    ----------
    axes : matplotlib Axes object
        Axes object that will get slider controls

    slider : `bool`
        ?

    Returns
    -------
    ? : ?
        ?

    .. todo::
        improve documentation
    """

    # If No axes specified use current axes.
    if not axes:
        axes = plt.gca()
    fig = axes.get_figure()

    # Split up the current axes so there is space for a start and a stop button
    divider = make_axes_locatable(axes)
    pad = 0.1 # Padding between axes
    pad_size = Size.Fraction(pad, Size.AxesX(axes))

    # Define size of useful axes cells, 50% each in x 20% for buttons in y.
    xsize = Size.Fraction((1.-2.*pad)/3., Size.AxesX(axes))
    ysize = Size.Fraction((1.-2.*pad)/15., Size.AxesY(axes))

    # Set up grid, 3x3 with cells for padding.
    divider.set_horizontal([xsize, pad_size, xsize, pad_size, xsize])
    if slider:
        divider.set_vertical([ysize, pad_size, ysize, pad_size, Size.AxesY(axes)])
        bny = 2
    else:
        divider.set_vertical([ysize, pad_size, Size.AxesY(axes)])
        bny = 0

    # Main figure spans all horiz and is in the top (2) in vert.
    axes.set_axes_locator(divider.new_locator(0, len(divider.get_vertical())-1,
                                              nx1=-1))

    # Add two axes for buttons and make them 50/50 spilt at the bottom.
    bax1 = fig.add_axes((0.,0.,1.,1.))
    locator = divider.new_locator(nx=0, ny=bny)
    bax1.set_axes_locator(locator)
    bax2 = fig.add_axes((0.,0.,0.8,1.))
    locator = divider.new_locator(nx=2, ny=bny)
    bax2.set_axes_locator(locator)
    bax3 = fig.add_axes((0.,0.,0.7,1.))
    locator = divider.new_locator(nx=4, ny=bny)
    bax3.set_axes_locator(locator)

    start = widgets.Button(bax1, "Start")
    stop = widgets.Button(bax2, "Stop")
    step = widgets.Button(bax3, "Step")
    # Make dummy reference to prevent garbage collection
    bax1._button = start
    bax2._button = stop
    bax3._button = step

    if slider:
        bax4 = fig.add_axes((0.,0.,0.6,1.))
        locator = divider.new_locator(nx=0, ny=0, nx1=-1)
        bax4.set_axes_locator(locator)
        sframe = widgets.Slider(bax4, 'Frame', 0, 10, valinit=0, valfmt = '%i')
        bax4._slider = sframe

        return axes, bax1, bax2, bax3, bax4
    return axes, bax1, bax2, bax3
