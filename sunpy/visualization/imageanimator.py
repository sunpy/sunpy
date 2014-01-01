# -*- coding: utf-8 -*-
import copy

import matplotlib.pyplot as plt
import matplotlib.widgets as widgets

from mpl_toolkits.axes_grid1 import make_axes_locatable
import mpl_toolkits.axes_grid1.axes_size as Size

__all__ = ['DataExplorer']

class SliderPB(widgets.Slider):
    __doc__= widgets.Slider.__doc__

    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, **kwargs):

        widgets.Slider.__init__(self, ax, label, valmin, valmax, valinit=valinit, valfmt=valfmt,
                 closedmin=closedmin, closedmax=closedmax, slidermin=slidermin,
                 slidermax=slidermax, dragging=dragging, **kwargs)

        self.changed_args = {}

    def set_val(self, val):
        xy = self.poly.xy
        xy[2] = val, 1
        xy[3] = val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson:
            return
        for cid, func in self.observers.iteritems():
            func(val, *self.changed_args[cid])

    def on_changed(self, func, *args):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.changed_args[cid] = args
        self.cnt += 1
        return cid

class ButtonPB(widgets.Button):
    def __init__(self, ax, label, image=None,
                 color='0.85', hovercolor='0.95'):

        widgets.Button.__init__(self, ax, label, image=image,
                 color=color, hovercolor=hovercolor)

        self.clicked_args = {}

    def on_clicked(self, func, *args):
        """
        When the button is clicked, call this *func* with event

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.clicked_args[cid] = args
        self.cnt += 1
        return cid

    def _release(self, event):
        if self.ignore(event):
            return
        if event.canvas.mouse_grabber != self.ax:
            return
        event.canvas.release_mouse(self.ax)
        if not self.eventson:
            return
        if event.inaxes != self.ax:
            return
        for cid, func in self.observers.iteritems():
            func(event, *self.clicked_args[cid])

class DataExplorer(object):
    """
    Create a matplotlib backend independant data explorer

    Parameters
    ----------
    data: ndarray
        The data to be visualised > 2D

    image_axes: list
        The two axes that make the image

    fig: mpl.figure
        Figure to use

    axis_range: list of lists
        List of [min, max] pairs for each axis

    interval: int
        Animation interval in ms

    Extra keywords are passed to imshow.
    """
    button_labels = ["Stop"]
    num_buttons = len(button_labels)

    def __init__(self, data, image_axes=[-2,-1], fig=None,
                 axis_range=None, interval=200, **kwargs):
        if not fig:
            fig = plt.figure()
        self.fig = fig

        self.data = data
        self.interval = interval
        self.naxis = data.ndim
        self.num_sliders = self.naxis - 2
        if len(image_axes) != 2:
            raise ValueError("There can only be two spatial axes")

        all_axes_f = range(self.naxis)
        self.image_axes = [all_axes_f[i] for i in image_axes]

        all_axes = range(self.naxis)
        [all_axes.remove(x) for x in self.spatial_axes]
        slider_axes = all_axes
        #Input Checking
        if len(slider_axes) != self.num_sliders:
            raise ValueError("Specific the same number of axes as sliders!")
        self.slider_axes = slider_axes

        ax = self.slider_axes + self.spatial_axes
        ax.sort()
        if ax != range(self.naxis):
            raise ValueError("spatial_axes and sider_axes mismatch")

        if not axis_range:
            axis_range = [[0, i] for i in self.data.shape]
        self.axis_range = axis_range

#==============================================================================
#         Begin Plotting etc.
#==============================================================================
        #Set a blank timer
        self.timer = None

        #Set up axes
        self._make_axes_grid()
        self._add_widgets()

        #create data slice
        self.spatial_slice = [slice(None)]*self.naxis
        for i in self.slider_axes:
            self.spatial_slice[i] = 0

        imshow_args = {'interpolation':'nearest',
                       'origin':'lower'}

        imshow_args.update(kwargs)
        self.im = self.axes.imshow(self.data[self.spatial_slice], **imshow_args)

        #Set the current axes to the main axes so commands like plt.ylabel() work.
        plt.sca(self.axes)

    def _get_slice(self, i, n):
        """
        Slice an array, i'th element on n'th axes

        Parameters
        ----------
        i: int
            The element to select
        n: int
            The axis along which to index the i'th element
        """
        nax = self.naxis
        arr_slice = [slice(None)]*nax
        arr_slice[n] = i
        return arr_slice

    def _updatefig(self, ax_slice):
        self.im.set_array(self.data[ax_slice])

    def _make_axes_grid(self):
        self.axes = self.fig.add_subplot(111)

        #Split up the current axes so there is space for a start and a stop button
        self.divider = make_axes_locatable(self.axes)
        pad = 0.01 # Padding between axes
        pad_size = Size.Fraction(pad, Size.AxesX(self.axes))
        large_pad_size = Size.Fraction(0.1, Size.AxesX(self.axes))

        #Define size of useful axes cells, 50% each in x 20% for buttons in y.
        small_x = Size.Fraction((1.-2.*pad)/10, Size.AxesX(self.axes))
        xsize = Size.Fraction((1.-2.*pad)/self.num_buttons, Size.AxesX(self.axes))
        ysize = Size.Fraction((1.-2.*pad)/15., Size.AxesY(self.axes))

        #Set up grid, 3x3 with cells for padding.
        self.divider.set_horizontal([xsize] + [pad_size, xsize]*(self.num_buttons-1)
                             + [Size.Fraction(0.1, Size.AxesY(self.axes)), small_x])

        vert = [ysize, pad_size] * self.num_sliders + [ysize, large_pad_size, Size.AxesY(self.axes)]
        self.divider.set_vertical(vert)
        self.button_ny = len(vert) - 3

        #Main figure spans all horiz and is in the top (2) in vert.
        self.axes.set_axes_locator(self.divider.new_locator(0, len(vert)-1,
                                                  nx1=-1))

    def _add_widgets(self):
        button_func = [self._stop_play]
        self.buttons = []
        for i in range(0,self.num_buttons):
            x = i*2
            #The i+1/10. is a bug that if you make two axes directly ontop of
            #one another then the divider doesn't work.
            self.buttons.append(self.fig.add_axes((0.,0.,0.+i/10.,1.)))
            locator = self.divider.new_locator(nx=x, ny=self.button_ny)
            self.buttons[-1].set_axes_locator(locator)
            self.buttons[-1]._button = widgets.Button(self.buttons[-1], self.button_labels[i])
            self.buttons[-1]._button.on_clicked(button_func[i])

        self.sliders = []
        self.radio = []
        for i in range(self.num_sliders):
            x = i * 2
            self.sliders.append(self.fig.add_axes((0.,0.,0.01+i/10.,1.)))
            locator = self.divider.new_locator(nx=0, ny=x, nx1=-3)
            self.sliders[-1].set_axes_locator(locator)
            sframe = SliderPB(self.sliders[-1], "%i"%i,
                                    self.axis_range[self.slider_axes[i]][0],
                                    self.axis_range[self.slider_axes[i]][1]-1,
                                    valinit=0, valfmt = '%i')
            sframe.on_changed(self._slider_changed, sframe)
            sframe.axes_num = self.slider_axes[i]
            sframe.cval = sframe.val
            self.sliders[-1]._slider = sframe

            self.radio.append(self.fig.add_axes((0., 0., 0.05+x/10., 1.)))
            locator = self.divider.new_locator(nx=2 + 2*(self.num_buttons-1), ny=x)
            self.radio[-1].set_axes_locator(locator)
            rdo = ButtonPB(self.radio[-1],'>>')
            rdo.on_clicked(self._start_play, sframe)
            self.radio[-1]._radio = rdo

    def _start_play(self, event, slider=None):
        if slider is None:
            slider = self.sliders[-1]._slider
        if not self.timer:
            self.timer = self.fig.canvas.new_timer()
            self.timer.interval = self.interval
            self.timer.add_callback(self._step, slider)
            self.timer.start()

    def _stop_play(self, event):
        if self.timer:
            self.timer.remove_callback(self._step)
            self.timer = None

    def _step(self, slider):
        s = slider
        if s.val >= s.valmax:
            s.set_val(s.valmin)
        else:
            s.set_val(s.val+1)

    def _slider_changed(self, val, slider):
        val = int(val)
        ax = slider.axes_num
        ax_slice = copy.copy(self.spatial_slice)
        ax_slice[ax] = val
        if val != slider.cval:
            self._updatefig(ax_slice)
            slider.cval = val


    def label_slider(self, i, label):
        """
        Change the Slider label

        Parameters
        ----------
        i: int
            The index of the slider to change (0 is bottom)
        label: str
            The label to set
        """
        self.sliders[i]._slider.label.set_text(label)