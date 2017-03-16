# -*- coding: utf-8 -*-

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import matplotlib.animation as mplanim

from mpl_toolkits.axes_grid1 import make_axes_locatable
import mpl_toolkits.axes_grid1.axes_size as Size

from sunpy.extern.six.moves import range

__all__ = ['BaseFuncAnimator', 'ImageAnimator']


class SliderPB(widgets.Slider):
    __doc__ = widgets.Slider.__doc__

    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, **kwargs):

        widgets.Slider.__init__(self, ax, label, valmin, valmax,
                                valinit=valinit, valfmt=valfmt,
                                closedmin=closedmin, closedmax=closedmax,
                                slidermin=slidermin, slidermax=slidermax,
                                dragging=dragging, **kwargs)
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
        for cid, func in self.observers.items():
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
        for cid, func in self.observers.items():
            func(event, *self.clicked_args[cid])


class BaseFuncAnimator(object):
    """
    Create a matplotlib backend independent data explorer which allows
    definition of figure update functions for each slider.

    The following keyboard shortcuts are defined in the viewer:

    - 'left': previous step on active slider
    - 'right': next step on active slider
    - 'top': change the active slider up one
    - 'bottom': change the active slider down one
    - 'p': play/pause active slider

    This viewer can have user defined buttons added by specifying the labels
    and functions called when those buttons are clicked as keyword arguments.

    To make this class useful the subclass must implement `_plot_start_image`
    which must define a `self.im` attribute which is an instance of AxesImage

    Parameters
    ----------
    data: iterable
        Some arbitrary data

    slider_functions: list
        A list of functions to call when that slider is changed.
        These functions will have `val`, the axes image object and the slider
        widget instance passed to them, i.e.: update_slider(val, im, slider)

    slider_ranges: list
        list of [min,max] pairs to set the ranges for each slider or an array
        of values for all points of the slider.
        (The slider update function decides which to support.)

    fig: mpl.figure
        Figure to use

    interval: int
        Animation interval in ms

    colorbar: bool
        Plot colorbar

    button_labels: list
        List of strings to label buttons

    button_func: list
        List of functions to map to the buttons

    Extra keywords are passed to imshow.
    """

    def __init__(self, data, slider_functions, slider_ranges, fig=None,
                 interval=200, colorbar=False, **kwargs):

        # Allow the user to specify the button func:
        self.button_func = kwargs.pop('button_func', [])
        self.button_labels = kwargs.pop('button_labels', [])
        self.num_buttons = len(self.button_labels)

        if not fig:
            fig = plt.figure()
        self.fig = fig

        self.data = data
        self.interval = interval
        self.if_colorbar = colorbar
        self.imshow_kwargs = kwargs

        if len(slider_functions) != len(slider_ranges):
            raise ValueError("You must specify the same number of functions "
                             "as extents")
        self.num_sliders = len(slider_functions)
        self.slider_functions = slider_functions
        self.slider_ranges = slider_ranges

        # Set active slider
        self.active_slider = 0

        # Set a blank timer
        self.timer = None

        # Set up axes
        self._make_axes_grid()
        self._add_widgets()
        self._set_active_slider(0)

        # Set the current axes to the main axes so commands like
        # plt.ylabel() work.
        plt.sca(self.axes)

        # Do Plot
        self.im = self.plot_start_image(self.axes)

        # Connect fig events
        self._connect_fig_events()

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

    def get_animation(self, axes=None, slider=0, startframe=0, endframe=None,
                      stepframe=1, **kwargs):
        """
        Return a matplotlib.animation.FuncAnimation instance for the selected
        slider.

        This will allow easy saving of the animation to a file.

        Parameters
        ----------

        slider: int
            The slider to animate along

        startframe: int
            The frame to start the animation

        endframe: int
            The frame to end the animation

        stepframe: int
            The step between frames
        """
        if not axes:
            axes = plt.gca()
        anim_fig = axes.get_figure()

        if endframe is None:
            endframe = self.slider_ranges[slider][1]

        im = self.plot_start_image(axes)

        anim_kwargs = {'frames': list(range(startframe, endframe, stepframe)),
                       'fargs': [im, self.sliders[slider]._slider]}
        anim_kwargs.update(kwargs)

        ani = mplanim.FuncAnimation(anim_fig, self.slider_functions[slider],
                                    **anim_kwargs)

        return ani

    def plot_start_image(self, ax):
        """
        This method creates the initial image on the mpl axes

        .. warning::
            This method needs to be implemented in subclasses

        Parameters
        ----------
        ax: mpl axes
            This is the axes on which to plot the image

        Returns
        -------
        AxesImage:
            An AxesImage object, the instance returned from a plt.imshow()
            command.
        """
        raise NotImplementedError("Please define your setup function")

    def _connect_fig_events(self):
        self.fig.canvas.mpl_connect('button_press_event', self._mouse_click)
        self.fig.canvas.mpl_connect('key_press_event', self._key_press)

    def _add_colorbar(self, im):
        self.colorbar = plt.colorbar(im, self.cax)

# =============================================================================
#   Figure event callback functions
# =============================================================================
    def _mouse_click(self, event):
        if event.inaxes in self.sliders:
            slider = self.sliders.index(event.inaxes)
            self._set_active_slider(slider)

    def _key_press(self, event):
        if event.key == 'left':
            self._previous(self.sliders[self.active_slider]._slider)
        elif event.key == 'right':
            self._step(self.sliders[self.active_slider]._slider)
        elif event.key == 'up':
            self._set_active_slider((self.active_slider+1) % self.num_sliders)
        elif event.key == 'down':
            self._set_active_slider((self.active_slider-1) % self.num_sliders)
        elif event.key == 'p':
            self._click_slider_button(event, self.slider_buttons[self.active_slider]._button,
                                      self.sliders[self.active_slider]._slider)

# =============================================================================
#   Active Slider methods
# =============================================================================
    def _set_active_slider(self, ind):
        self._dehighlight_slider(self.active_slider)
        self._highliget_slider(ind)
        self.active_slider = ind

    def _highliget_slider(self, ind):
        ax = self.sliders[ind]
        [a.set_linewidth(2.0) for n, a in ax.spines.items()]
        self.fig.canvas.draw()

    def _dehighlight_slider(self, ind):
        ax = self.sliders[ind]
        [a.set_linewidth(1.0) for n, a in ax.spines.items()]
        self.fig.canvas.draw()

# =============================================================================
#   Build the figure and place the widgets
# =============================================================================
    def _get_main_axes(self):
        """ Allow replacement of main axes by subclassing """
        return self.fig.add_subplot(111)

    def _make_axes_grid(self):
        self.axes = self._get_main_axes()

        # Split up the current axes so there is space for start & stop buttons
        self.divider = make_axes_locatable(self.axes)
        pad = 0.01  # Padding between axes
        pad_size = Size.Fraction(pad, Size.AxesX(self.axes))
        large_pad_size = Size.Fraction(0.1, Size.AxesY(self.axes))

        # Define size of useful axes cells, 50% each in x 20% for buttons in y.
        small_x = Size.Fraction((1.-2.*pad)/10, Size.AxesX(self.axes))
        ysize = Size.Fraction((1.-2.*pad)/15., Size.AxesY(self.axes))

        # Set up grid, 3x3 with cells for padding.
        if self.num_buttons > 0:
            xsize = Size.Fraction((1.-2.*pad)/self.num_buttons, Size.AxesX(self.axes))
            horiz = [xsize] + [pad_size, xsize]*(self.num_buttons-1) + \
                    [Size.Fraction(0.1, Size.AxesY(self.axes)), small_x]
            vert = [ysize, pad_size] * self.num_sliders + \
                   [large_pad_size, large_pad_size, Size.AxesY(self.axes)]
        else:
            vert = [ysize, pad_size] * self.num_sliders + \
                   [large_pad_size, Size.AxesY(self.axes)]
            horiz = [Size.Fraction(0.8, Size.AxesX(self.axes))] + \
                    [Size.Fraction(0.1, Size.AxesX(self.axes))]*2

        self.divider.set_horizontal(horiz)
        self.divider.set_vertical(vert)
        self.button_ny = len(vert) - 3

        # If we are going to add a colorbar it'll need an axis next to the plot
        if self.if_colorbar:
            nx1 = -3
            self.cax = self.fig.add_axes((0., 0., 0.141, 1.))
            locator = self.divider.new_locator(nx=-2, ny=len(vert)-1, nx1=-1)
            self.cax.set_axes_locator(locator)
        else:
            # Main figure spans all horiz and is in the top (2) in vert.
            nx1 = -1

        self.axes.set_axes_locator(
            self.divider.new_locator(nx=0, ny=len(vert)-1, nx1=nx1))

    def _add_widgets(self):
        self.buttons = []
        for i in range(0, self.num_buttons):
            x = i*2
            # The i+1/10. is a bug that if you make two axes directly on top of
            # one another then the divider doesn't work.
            self.buttons.append(self.fig.add_axes((0., 0., 0.+i/10., 1.)))
            locator = self.divider.new_locator(nx=x, ny=self.button_ny)
            self.buttons[-1].set_axes_locator(locator)
            self.buttons[-1]._button = widgets.Button(self.buttons[-1],
                                                      self.button_labels[i])
            self.buttons[-1]._button.on_clicked(self.button_func[i])

        self.sliders = []
        self.slider_buttons = []
        for i in range(self.num_sliders):
            x = i * 2
            self.sliders.append(self.fig.add_axes((0., 0., 0.01+i/10., 1.)))
            if self.num_buttons == 0:
                nx1 = 1
            else:
                nx1 = -3
            locator = self.divider.new_locator(nx=0, ny=x, nx1=nx1)
            self.sliders[-1].set_axes_locator(locator)
            sframe = SliderPB(self.sliders[-1], "{slide:d}".format(slide=i),
                              self.slider_ranges[i][0],
                              self.slider_ranges[i][-1]-1,
                              valinit=self.slider_ranges[i][0],
                              valfmt='%4.1f')
            sframe.on_changed(self._slider_changed, sframe)
            sframe.slider_ind = i
            sframe.cval = sframe.val
            self.sliders[-1]._slider = sframe

            self.slider_buttons.append(
                self.fig.add_axes((0., 0., 0.05+x/10., 1.)))
            if self.num_buttons == 0:
                nx = 2
            else:
                nx = 2 + 2*(self.num_buttons-1)
            locator = self.divider.new_locator(nx=nx, ny=x)

            self.slider_buttons[-1].set_axes_locator(locator)
            butt = ButtonPB(self.slider_buttons[-1], ">")
            butt.on_clicked(self._click_slider_button, butt, sframe)
            butt.clicked = False
            self.slider_buttons[-1]._button = butt

# =============================================================================
#   Widget callbacks
# =============================================================================
    def _click_slider_button(self, event, button, slider):
        self._set_active_slider(slider.slider_ind)
        if button.clicked:
            self._stop_play(event)
            button.clicked = False
            button.label.set_text(">")
        else:
            button.clicked = True
            self._start_play(event, button, slider)
            button.label.set_text("||")
        self.fig.canvas.draw()

    def _start_play(self, event, button, slider):
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
        self.fig.canvas.draw()

    def _previous(self, slider):
        s = slider
        if s.val <= s.valmin:
            s.set_val(s.valmax)
        else:
            s.set_val(s.val-1)
        self.fig.canvas.draw()

    def _slider_changed(self, val, slider):
        self.slider_functions[slider.slider_ind](val, self.im, slider)


class ArrayAnimator(BaseFuncAnimator):
    """
    Create a matplotlib backend independent data explorer

    The following keyboard shortcuts are defined in the viewer:

    - 'left': previous step on active slider
    - 'right': next step on active slider
    - 'top': change the active slider up one
    - 'bottom': change the active slider down one
    - 'p': play/pause active slider

    This viewer can have user defined buttons added by specifying the labels
    and functions called when those buttons are clicked as keyword arguments.

    Parameters
    ----------
    data: ndarray
        The data to be visualised

    image_axes: list
        The axes that make the image

    fig: mpl.figure
        Figure to use

    axis_range: list of physical coordinates for array or None
        If None array indices will be used for all axes.
        If a list it should contain one element for each axis of the numpy array.
        For the image axes a [min, max] pair should be specified which will be
        passed to :func:`matplotlib.pyplot.imshow` as extent.
        For the slider axes a [min, max] pair can be specified or an array the
        same length as the axis which will provide all values for that slider.
        If None is specified for an axis then the array indices will be used
        for that axis.

    interval: int
        Animation interval in ms

    colorbar: bool
        Plot colorbar

    button_labels: list
        List of strings to label buttons

    button_func: list
        List of functions to map to the buttons

    """

    def __init__(self, data, image_axes=[-2, -1], axis_range=None, **kwargs):

        all_axes = list(range(self.naxis))
        # Handle negative indexes
        self.image_axes = [all_axes[i] for i in image_axes]

        slider_axes = list(range(self.naxis))
        for x in self.image_axes:
            slider_axes.remove(x)

        if len(slider_axes) != self.num_sliders:
            raise ValueError("Specific the same number of axes as sliders!")
        self.slider_axes = slider_axes

        # Verify that combined slider_axes and image_axes make all axes
        ax = self.slider_axes + self.image_axes
        ax.sort()
        if ax != list(range(self.naxis)):
            raise ValueError("spatial_axes and slider_axes mismatch")

        self.axis_range = self._sanitize_axis_range(axis_range, data)

        # create data slice
        self.frame_slice = [slice(None)] * self.naxis
        for i in self.slider_axes:
            self.frame_slice[i] = 0

        base_kwargs = {'slider_functions': [self.update_plot] * self.num_sliders,
                       'slider_ranges': [self.axis_range[i] for i in self.slider_axes]}
        base_kwargs.update(kwargs)
        BaseFuncAnimator.__init__(self, data, **base_kwargs)

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

    def _sanitize_axis_range(self, axis_range, data):
        """
        This method takes the various allowed values of axis_range and returns
        them in a standardized way for the rest of the class to use.

        The outputted axis range describes the physical coordinates of the
        array axes.

        The allowed values of axis range is either None or a list.
        If axis_range is None then all axis are assumed to be not scaled and
        use array indices.

        Where axis_range is a list it must have the same length as the number
        of axis as the array and each element must be one of the following:

            * None: Build a min,max pair or linspace array of array indices
            * [min, max]: leave for image axes or convert to a array for slider axes
            (from min to max in axis length steps)
            * [min, max] pair where min == max: convert to array indies min,max pair or array.
            * array of axis length, check that it was passed for a slider axes and do nothing
            if it was, error if it is not.
        """
        # If no axis range at all make it all [min,max] pairs
        if axis_range is None:
            axis_range = [[0, i] for i in data.shape]

        # need the same number of axis ranges as axes
        if len(axis_range) != data.ndim:
            raise ValueError("Length of axis_range must equal number of axes")

        # For each axis validate and translate the axis_range
        for i, d in enumerate(data.shape):
            # If [min,max] pair or None
            if axis_range[i] is None or len(axis_range[i]) == 2:
                # If min==max or None
                if axis_range[i] is None or axis_range[i][0] == axis_range[i][1]:
                    if i in self.slider_axes:
                        axis_range[i] = np.linspace(0, d, d)
                    else:
                        axis_range[i] = [0, d]
                        # min max pair for slider axes should be converted
                        # to an array
                elif i in self.slider_axes:
                    axis_range[i] = np.linspace(axis_range[i][0], axis_range[i][1], d)

            # If we have a whole list of values for the axis, make sure we are a slider axis.
            elif len(axis_range[i]) == d:
                axis_range[i] = np.array(axis_range[i])

            # panic
            else:
                raise ValueError("axis_range must be None or a list with length equal to number "
                                 "of axes in data whose elements are either None, [min,max], "
                                 "or a list/array of same length as the plot/image axis of data.")
        return axis_range


class ImageAnimator(ArrayAnimator):
    """Create a matplotlib backend independent data explorer for 2D images.

    The following keyboard shortcuts are defined in the viewer:

    - 'left': previous step on active slider
    - 'right': next step on active slider
    - 'top': change the active slider up one
    - 'bottom': change the active slider down one
    - 'p': play/pause active slider

    This viewer can have user defined buttons added by specifying the labels
    and functions called when those buttons are clicked as keyword arguments.

    Parameters
    ----------
    data: ndarray
        The data to be visualised >= 2D

    image_axes: list
        The two axes that make the image

    fig: mpl.figure
        Figure to use

    axis_range: list of physical coordinates for array or None
        If None array indices will be used for all axes.
        If a list it should contain one element for each axis of the numpy array.
        For the image axes a [min, max] pair should be specified which will be
        passed to :func:`matplotlib.pyplot.imshow` as extent.
        For the slider axes a [min, max] pair can be specified or an array the
        same length as the axis which will provide all values for that slider.
        If None is specified for an axis then the array indices will be used
        for that axis.

    interval: int
        Animation interval in ms

    colorbar: bool
        Plot colorbar

    button_labels: list
        List of strings to label buttons

    button_func: list
        List of functions to map to the buttons

    Extra keywords are passed to imshow.

    """

    def __init__(self, data, image_axes=[-2, -1], axis_range=None, **kwargs):
        # Check that number of axes is 2.
        if len(image_axes) != 2:
            raise ValueError("There can only be two spatial axes")
        # Define number of slider axes.
        self.naxis = data.ndim
        self.num_sliders = self.naxis-2
        # Run init for parent class
        super(ImageAnimator, self).__init__(data, image_axes=image_axes,
                                            axis_range=axis_range, **kwargs)

    def plot_start_image(self, ax):
        # Create extent arg
        extent = []
        # reverse because numpy is in y-x and extent is x-y
        for i in self.image_axes[::-1]:
            extent += self.axis_range[i]

        imshow_args = {'interpolation': 'nearest',
                       'origin': 'lower',
                       'extent': extent,
                       }

        imshow_args.update(self.imshow_kwargs)
        im = ax.imshow(self.data[self.frame_slice], **imshow_args)
        if self.if_colorbar:
            self._add_colorbar(im)

        return im

    def update_plot(self, val, im, slider):
        val = int(val)
        ax_ind = self.slider_axes[slider.slider_ind]
        ind = np.argmin(np.abs(self.axis_range[ax_ind] - val))
        self.frame_slice[ax_ind] = ind
        if val != slider.cval:
            im.set_array(self.data[self.frame_slice])
            slider.cval = val

    def _sanitize_axis_range(self, axis_range, data):
        """
        This method takes the various allowed values of axis_range and returns
        them in a standardized way for the rest of the class to use.

        The outputted axis range describes the physical coordinates of the
        array axes.

        The allowed values of axis range is either None or a list.
        If axis_range is None then all axis are assumed to be not scaled and
        use array indices.

        Where axis_range is a list it must have the same length as the number
        of axis as the array and each element must be one of the following:

            * None: Build a min,max pair or linspace array of array indices
            * [min, max]: leave for image axes or convert to a array for slider axes
            (from min to max in axis length steps)
            * [min, max] pair where min == max: convert to array indies min,max pair or array.
            * array of axis length, check that it was passed for a slider axes and do nothing
            if it was, error if it is not.
        """
        # Run super class's version of this function.
        axis_range = super(ImageAnimator, self)._sanitize_axis_range(axis_range, data)
        # Check that axis ranges which are array type are only for slider axes.
        for i, d in enumerate(data.shape):
            if len(axis_range[i]) == d:
                if i not in self.slider_axes:
                    raise ValueError("Slider axes mis-match, non-slider axes need [min,max] pairs")
        return axis_range



class LineAnimator(ArrayAnimator):
    """Create a matplotlib backend independent data explorer for 1D plots.

    The following keyboard shortcuts are defined in the viewer:

    - 'left': previous step on active slider
    - 'right': next step on active slider
    - 'top': change the active slider up one
    - 'bottom': change the active slider down one
    - 'p': play/pause active slider

    This viewer can have user defined buttons added by specifying the labels
    and functions called when those buttons are clicked as keyword arguments.

    Parameters
    ----------
    data: ndarray
        The y-axis data to be visualised = 2D

    data_yaxis: `int`
        The axis used to plot against xdata.

    fig: `matplotlib.figure`
        Figure to use

    axis_range: `list` of physical coordinates for array or None
        If None array indices will be used for all axes.
        X-axis values must be supplied (if desired) as an array in the
        element of axis_range corresponding to the data_yaxis in the
        data input arg, i.e. xdata = axis_range[data_yaxis].
        Also, the number of x-axis values must correspond to the number
        of y-axis values, i.e. len(axis_range[data_yaxis]) must equal
        len(data[data_yaxis]).
        For the slider axes a [min, max] pair can be specified or an array the
        same length as the axis which will provide all values for that slider.
        If None is specified for an axis then the array indices will be used
        for that axis.

    xlabel: `str`
        Label of x-axis of plot.

    ylabel: `str`
        Label of y-axis of plot.

    xlim: `tuple`
        Limits of x-axis of plot.

    ylim: `tuple`
        Limits of y-axis of plot.

    interval: `int`
        Animation interval in ms

    button_labels: `list`
        List of strings to label buttons

    button_func: `list`
        List of functions to map to the buttons

    Extra keywords are passed to plot.

    """

    def __init__(self, data, data_yaxis=0, axis_range=None, ylabel=None, xlabel=None,
                 xlim=None, ylim=None, **kwargs):
        # Check inputs.
        data_yaxis = int(data_yaxis)
        if data_yaxis != 0 and data_yaxis != 1:
            raise ValueError("data_yaxis must be either 0 or 1 referring to "
                             "the axis of data to be used for a single plot.")
        if data.ndim != 2:
            raise ValueError("data must have two dimensions.  One for data "
                             "for each single plot and one for time/iteration.")
        # Attach data to class.
        axis_range[data_yaxis] = np.array(axis_range[data_yaxis])
        self.xdata = axis_range[data_yaxis]
        if ylim is None:
            ylim = (data.min(), data.max())
        if xlim is None:
            xlim = (self.xdata.min(), self.xdata.max())
        self.ylim = ylim
        self.xlim = xlim
        self.xlabel = xlabel
        self.ylabel = ylabel
        # Define number of slider axes.
        self.naxis = data.ndim
        self.num_sliders = self.naxis-1
        # Run init for base class
        super(LineAnimator, self).__init__(data, image_axes=[data_yaxis],
                                           axis_range=axis_range, **kwargs)

    def plot_start_image(self, ax):
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        if self.xlabel is not None:
            ax.set_xlabel(self.xlabel)
        if self.ylabel is not None:
            ax.set_ylabel(self.ylabel)
        plot_args = {}
        plot_args.update(self.imshow_kwargs)
        line, = ax.plot(self.xdata, self.data[self.frame_slice], **plot_args)
        return line

    def update_plot(self, val, line, slider):
        val = int(val)
        ax_ind = self.slider_axes[slider.slider_ind]
        ind = np.argmin(np.abs(self.axis_range[ax_ind] - val))
        self.frame_slice[ax_ind] = ind
        if val != slider.cval:
            line.set_ydata(self.data[self.frame_slice])
            slider.cval = val
