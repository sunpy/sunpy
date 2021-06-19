import abc
from functools import partial

import matplotlib.animation as mplanim
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import mpl_toolkits.axes_grid1.axes_size as Size
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

import astropy.units as u

__all__ = ['BaseFuncAnimator', 'ArrayAnimator']


class BaseFuncAnimator:
    """
    Create a Matplotlib backend independent data explorer which allows
    definition of figure update functions for each slider.

    The following keyboard shortcuts are defined in the viewer:

    * 'left': previous step on active slider.
    * 'right': next step on active slider.
    * 'top': change the active slider up one.
    * 'bottom': change the active slider down one.
    * 'p': play/pause active slider.

    This viewer can have user defined buttons added by specifying the labels
    and functions called when those buttons are clicked as keyword arguments.

    To make this class useful the subclass must implement ``_plot_start_image``
    which must define a ``self.im`` attribute which is an instance of
    `matplotlib.image.AxesImage`.

    Parameters
    ----------
    data: `iterable`
        Some arbitrary data.
    slider_functions: `list`
        A list of functions to call when that slider is changed.
        These functions will have ``val``, the axes image object and the slider
        widget instance passed to them, e.g., ``update_slider(val, im, slider)``
    slider_ranges: `list`
        A list of ``[min,max]`` pairs to set the ranges for each slider or an array
        of values for all points of the slider.
        (The slider update function decides which to support.)
    fig: `matplotlib.figure.Figure`, optional
        `~matplotlib.figure.Figure` to use. Defaults to `None`, in which case a new
        figure is created.
    interval: `int`, optional
        Animation interval in milliseconds. Defaults to 200.
    colorbar: `bool`, optional
        Plot a colorbar. Defaults to `False`.
    button_labels: `list`, optional
        A list of strings to label buttons. Defaults to `None`. If `None`
        and ``button_func`` is specified, it will default to the
        names of the functions.
    button_func: `list`, optional
        A list of functions to map to the buttons. These functions are called
        with two arguments, ``(animator, event)`` where the first argument is
        the animator object, and the second is a
        `matplotlib.backend_bases.MouseEvent` object.
        Defaults to `None`.
    slider_labels: `list`, optional
        A list of labels to draw in the slider, must be the same length as
        ``slider_functions``.

    Notes
    -----
    Extra keywords are passed to `matplotlib.pyplot.imshow`.
    """

    def __init__(self, data, slider_functions, slider_ranges, fig=None,
                 interval=200, colorbar=False, button_func=None, button_labels=None,
                 start_image_func=None, slider_labels=None, **kwargs):

        # Allow the user to specify the button func:
        self.button_func = button_func or []
        if button_func and not button_labels:
            button_labels = [a.__name__ for a in button_func]
        self.button_labels = button_labels or []
        self.num_buttons = len(self.button_func)

        if not fig:
            fig = plt.figure()
        self.fig = fig

        self.data = data
        self.interval = interval
        self.if_colorbar = colorbar
        self.imshow_kwargs = kwargs

        if len(slider_functions) != len(slider_ranges):
            raise ValueError("slider_functions and slider_ranges must be the same length.")

        if slider_labels is not None:
            if len(slider_labels) != len(slider_functions):
                raise ValueError("slider_functions and slider_labels must be the same length.")

        self.num_sliders = len(slider_functions)
        self.slider_functions = slider_functions
        self.slider_ranges = slider_ranges
        self.slider_labels = slider_labels or [''] * len(slider_functions)

        # Set active slider
        self.active_slider = 0

        # Set a blank timer
        self.timer = None

        # Set up axes
        self.axes = None
        self._make_axes_grid()
        self._add_widgets()
        self._set_active_slider(0)

        # Set the current axes to the main axes so commands like plt.ylabel() work.
        #
        # Only do this if figure has a manager, so directly constructed figures
        # (ie. via matplotlib.figure.Figure()) work.
        if hasattr(self.fig.canvas, "manager") and self.fig.canvas.manager is not None:
            plt.sca(self.axes)

        # Do Plot
        self.im = self.plot_start_image(self.axes)

        # Connect fig events
        self._connect_fig_events()

    def label_slider(self, i, label):
        """
        Change the slider label.

        Parameters
        ----------
        i: `int`
            The index of the slider to change (0 is bottom).
        label: `str`
            The label to set.
        """
        self.sliders[i]._slider.label.set_text(label)

    def get_animation(self, axes=None, slider=0, startframe=0, endframe=None,
                      stepframe=1, **kwargs):
        """
        Return a `~matplotlib.animation.FuncAnimation` instance for the
        selected slider.

        This will allow easy saving of the animation to a file.

        Parameters
        ----------
        axes: `matplotlib.axes.Axes`, optional
            The `matplotlib.axes.Axes` to animate. Defaults to `None`, in which
            case the Axes associated with this animator are used. Passing a
            custom Axes can be useful if you want to create the animation on
            a custom figure that is not the figure set up by this Animator.
        slider: `int`, optional
            The slider to animate along. Defaults to 0.
        startframe: `int`, optional
            The frame to start the animation. Defaults to 0.
        endframe: `int`, optional
            The frame to end the animation. Defaults to `None`.
        stepframe: `int`, optional
            The step between frames. Defaults to 1.

        Notes
        -----
        Extra keywords are passed to `matplotlib.animation.FuncAnimation`.
        """
        if not axes:
            axes = self.axes
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
        This method creates the initial image on the `matplotlib.axes.Axes`.

        .. warning::
            This method needs to be implemented in subclasses.

        Parameters
        ----------
        ax: `matplotlib.axes.Axes`
            This is the axes on which to plot the image.

        Returns
        -------
        `matplotlib.artist.Artist`
            The matplotlib object to be animated, this is usually either a
            `~matplotlib.image.AxesImage` object, or a
            `~matplotlib.lines.Line2D`.
        """
        raise NotImplementedError("Please define this function.")

    def _connect_fig_events(self):
        self.fig.canvas.mpl_connect('button_press_event', self._mouse_click)
        self.fig.canvas.mpl_connect('key_press_event', self._key_press)

    def _add_colorbar(self, im):
        self.colorbar = self.fig.colorbar(im, self.cax)

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
        self._highlight_slider(ind)
        self.active_slider = ind

    def _highlight_slider(self, ind):
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
    def _setup_main_axes(self):
        """
        Allow replacement of main axes by subclassing.
        This method must set the ``axes`` attribute.
        """
        if self.axes is None:
            self.axes = self.fig.add_subplot(111)

    def _make_axes_grid(self):
        self._setup_main_axes()

        # Split up the current axes so there is space for start & stop buttons
        self.divider = make_axes_locatable(self.axes)
        pad = 0.01  # Padding between axes
        pad_size = Size.Fraction(pad, Size.AxesX(self.axes))
        large_pad_size = Size.Fraction(0.1, Size.AxesY(self.axes))

        button_grid = max((7, self.num_buttons))

        # Define size of useful axes cells, 50% each in x 20% for buttons in y.
        ysize = Size.Fraction((1.-2.*pad)/15., Size.AxesY(self.axes))
        xsize = Size.Fraction((1.-2.*pad)/button_grid, Size.AxesX(self.axes))

        # Set up grid, 3x3 with cells for padding.
        if self.num_buttons > 0:
            horiz = [xsize] + [pad_size, xsize]*(button_grid-1)
            vert = [ysize, pad_size] * self.num_sliders + \
                   [large_pad_size, large_pad_size, Size.AxesY(self.axes)]
        else:
            vert = [ysize, large_pad_size] * self.num_sliders + \
                   [large_pad_size, Size.AxesY(self.axes)]
            horiz = [Size.Fraction(0.1, Size.AxesX(self.axes))] + \
                    [Size.Fraction(0.05, Size.AxesX(self.axes))] + \
                    [Size.Fraction(0.65, Size.AxesX(self.axes))] + \
                    [Size.Fraction(0.1, Size.AxesX(self.axes))] + \
                    [Size.Fraction(0.1, Size.AxesX(self.axes))]

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
            x = i * 2
            # The i+1/10. is a bug that if you make two axes directly on top of
            # one another then the divider doesn't work.
            self.buttons.append(self.fig.add_axes((0., 0., 0.+i/10., 1.)))
            locator = self.divider.new_locator(nx=x, ny=self.button_ny)
            self.buttons[-1].set_axes_locator(locator)
            self.buttons[-1]._button = widgets.Button(self.buttons[-1],
                                                      self.button_labels[i])
            self.buttons[-1]._button.on_clicked(partial(self.button_func[i], self))

        self.sliders = []
        self.slider_buttons = []
        for i in range(self.num_sliders):
            y = i * 2
            self.sliders.append(self.fig.add_axes((0., 0., 0.01+i/10., 1.)))
            if self.num_buttons == 0:
                nx1 = 3
            else:
                nx1 = -2
            locator = self.divider.new_locator(nx=2, ny=y, nx1=nx1)
            self.sliders[-1].set_axes_locator(locator)
            self.sliders[-1].text(0.5, 0.5, self.slider_labels[i],
                                  transform=self.sliders[-1].transAxes,
                                  horizontalalignment="center",
                                  verticalalignment="center")

            sframe = widgets.Slider(self.sliders[-1], "",
                                    self.slider_ranges[i][0],
                                    self.slider_ranges[i][-1]-1,
                                    valinit=self.slider_ranges[i][0],
                                    valfmt='%4.1f')
            sframe.on_changed(partial(self._slider_changed, slider=sframe))
            sframe.slider_ind = i
            sframe.cval = sframe.val
            self.sliders[-1]._slider = sframe

            self.slider_buttons.append(
                self.fig.add_axes((0., 0., 0.05+y/10., 1.)))
            locator = self.divider.new_locator(nx=0, ny=y)

            self.slider_buttons[-1].set_axes_locator(locator)
            butt = widgets.Button(self.slider_buttons[-1], ">")
            butt.on_clicked(partial(self._click_slider_button, button=butt, slider=sframe))
            butt.clicked = False
            self.slider_buttons[-1]._button = butt

# =============================================================================
#   Widget callbacks
# =============================================================================
    def _slider_changed(self, val, slider):
        self.slider_functions[slider.slider_ind](val, self.im, slider)

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


class ArrayAnimator(BaseFuncAnimator, metaclass=abc.ABCMeta):
    """
    Create a Matplotlib backend independent data explorer.

    The following keyboard shortcuts are defined in the viewer:

    * 'left': previous step on active slider.
    * 'right': next step on active slider.
    * 'top': change the active slider up one.
    * 'bottom': change the active slider down one.
    * 'p': play/pause active slider.

    This viewer can have user defined buttons added by specifying the labels
    and functions called when those buttons are clicked as keyword arguments.

    Parameters
    ----------
    data: `numpy.ndarray`
        The data to be visualized.
    image_axes: `list`, optional
        A list of the axes order that make up the image.
    axis_ranges: `list` of physical coordinates for the `numpy.ndarray`, optional
        Defaults to `None` and array indices will be used for all axes.
        The `list` should contain one element for each axis of the `numpy.ndarray`.
        For the image axes a ``[min, max]`` pair should be specified which will be
        passed to `matplotlib.pyplot.imshow` as an extent.
        For the slider axes a ``[min, max]`` pair can be specified or an array the
        same length as the axis which will provide all values for that slider.

    Notes
    -----
    Extra keywords are passed to `~sunpy.visualization.animator.BaseFuncAnimator`.
    """

    def __init__(self, data, image_axes=[-2, -1], axis_ranges=None, **kwargs):

        all_axes = list(range(self.naxis))
        # Handle negative indexes
        self.image_axes = [all_axes[i] for i in image_axes]

        slider_axes = list(range(self.naxis))
        for x in self.image_axes:
            slider_axes.remove(x)

        if len(slider_axes) != self.num_sliders:
            raise ValueError("Number of sliders doesn't match the number of slider axes.")
        self.slider_axes = slider_axes

        # Verify that combined slider_axes and image_axes make all axes
        ax = self.slider_axes + self.image_axes
        ax.sort()
        if ax != list(range(self.naxis)):
            raise ValueError("Number of image and slider axes do not match total number of axes.")

        self.axis_ranges, self.extent = self._sanitize_axis_ranges(axis_ranges, data.shape)

        # create data slice
        self.frame_slice = [slice(None)] * self.naxis
        for i in self.slider_axes:
            self.frame_slice[i] = 0

        slider_functions = kwargs.pop("slider_functions", [])
        slider_ranges = kwargs.pop("slider_ranges", [])
        base_kwargs = {
            'slider_functions': ([self.update_plot] * self.num_sliders) + slider_functions,
            'slider_ranges': [[0, dim] for dim in np.array(data.shape)[self.slider_axes]] + slider_ranges
        }
        self.num_sliders = len(base_kwargs["slider_functions"])
        base_kwargs.update(kwargs)
        super().__init__(data, **base_kwargs)

    @property
    def frame_index(self):
        """
        A tuple version of ``frame_slice`` to be used when indexing arrays.
        """
        return tuple(self.frame_slice)

    def label_slider(self, i, label):
        """
        Change the Slider label.

        Parameters
        ----------
        i: `int`
            The index of the slider to change (0 is bottom).
        label: `str`
            The label to set.
        """
        self.sliders[i]._slider.label.set_text(label)

    def _sanitize_axis_ranges(self, axis_ranges, data_shape):
        """
        This method takes the various allowed values of ``axis_ranges`` and
        returns them in a standardized way for the rest of the class to use.

        The outputted axis range describes the physical coordinates of the
        array axes.

        The allowed values of axis range is either `None` or a `list`.
        If ``axis_ranges`` is `None` then all axis are assumed to be not scaled and
        will use array indices.

        Where ``axis_ranges`` is a `list` it must have the same length as the number
        of axis as the array and each element must be one of the following:

        * `None`: Build a "min,max" pair or `numpy.linspace` array of array indices.
        * ``[min, max]``: Either leave for the image axes or convert to a array for slider axes
          (from min to max in axis length steps)
        * ``[min, max]`` pair where ``min == max``: convert to array indies "min,max" pair or array.
        * array of axis length, check that it was passed for a slider axes and do nothing
          if it was, error if it is not.
        * For slider axes: a function which maps from pixel to world value.
        """
        ndim = len(data_shape)
        # If no axis range at all make it all [min,max] pairs
        if axis_ranges is None:
            axis_ranges = [None] * ndim

        # need the same number of axis ranges as axes
        if len(axis_ranges) != ndim:
            raise ValueError("Length of axis_ranges must equal number of axes")

        # Define error message for incompatible axis_range input.
        def incompatible_axis_ranges_error_message(j): return \
            (f"Unrecognized format for {j}th entry in axis_ranges: {axis_ranges[j]}"
             "axis_ranges must be None, a ``[min, max]`` pair, or "
             "an array-like giving the edge values of each pixel, "
             "i.e. length must be length of axis + 1.")

        # If axis range not given, define a function such that the range goes
        # from -0.5 to number of pixels-0.5.  Thus, the center of the pixels
        # along the axis will correspond to integer values.
        def none_image_axis_range(j): return [-0.5, data_shape[j]-0.5]

        # For each axis validate and translate the axis_ranges. For image axes,
        # also determine the plot extent.  To do this, iterate through image and slider
        # axes separately.  Iterate through image axes in reverse order
        # because numpy is in y-x and extent is x-y.
        extent = []
        for i in self.image_axes[::-1]:
            if axis_ranges[i] is None:
                extent = extent + none_image_axis_range(i)
                axis_ranges[i] = np.array(none_image_axis_range(i))
            else:
                # Depending on length of axis_ranges[i], leave unchanged,
                # convert to pixel centers or raise an error due to incompatible format.
                axis_ranges[i] = np.asarray(axis_ranges[i])
                if len(axis_ranges[i]) == 2:
                    # Set extent.
                    extent += [axis_ranges[i][0], axis_ranges[i][-1]]
                elif axis_ranges[i].ndim == 1 and len(axis_ranges[i]) == data_shape[i]+1:
                    # If array of individual pixel edges supplied, first set extent
                    # from first and last pixel edge, then convert axis_ranges to pixel centers.
                    # The reason that pixel edges are required as input rather than centers
                    # is so that the plot extent can be derived from axis_ranges (above)
                    # and APIs using both [min, max] pair and manual definition of each pixel
                    # values can be unambiguously and simultanously supported.
                    extent += [axis_ranges[i][0], axis_ranges[i][-1]]
                    axis_ranges[i] = edges_to_centers_nd(axis_ranges[i], 0)
                elif axis_ranges[i].ndim == ndim and axis_ranges[i].shape[i] == data_shape[i]+1:
                    extent += [axis_ranges[i].min(), axis_ranges[i].max()]
                    axis_ranges[i] = edges_to_centers_nd(axis_ranges[i], i)
                else:
                    raise ValueError(incompatible_axis_ranges_error_message(i))

        # For each slider axis validate and translate the axis_ranges.
        def get_pixel_to_world_callable(array):
            def pixel_to_world(pixel):
                return array[pixel]
            return pixel_to_world

        for sidx in self.slider_axes:
            if axis_ranges[sidx] is None:
                # If axis range not supplied, set pixel center values as integers starting at 0.
                axis_ranges[sidx] = get_pixel_to_world_callable(np.arange(data_shape[sidx]))
            elif not callable(axis_ranges[sidx]):
                axis_ranges[sidx] = np.array(axis_ranges[sidx])
                if len(axis_ranges[sidx]) == 2:
                    # If axis range given as a min, max pair, derive the center of each pixel
                    # assuming they are equally spaced.

                    axis_ranges[sidx] = np.linspace(axis_ranges[sidx][0], axis_ranges[sidx][-1],
                                                    data_shape[sidx]+1)
                    axis_ranges[sidx] = get_pixel_to_world_callable(
                        edges_to_centers_nd(axis_ranges[sidx], sidx))
                elif axis_ranges[sidx].ndim == 1 and len(axis_ranges[sidx]) == data_shape[sidx]+1:
                    # If axis range given as 1D array of pixel edges (i.e. axis is independent),
                    # derive pixel centers.

                    axis_ranges[sidx] = get_pixel_to_world_callable(
                        edges_to_centers_nd(np.asarray(axis_ranges[sidx]), 0))
                elif axis_ranges[sidx].ndim == ndim and axis_ranges[sidx].shape[sidx] == data_shape[sidx]+1:
                    # If axis range given as array of pixel edges the same shape as
                    # the data array (i.e. axis is not independent), derive pixel centers.

                    axis_ranges[sidx] = get_pixel_to_world_callable(
                        edges_to_centers_nd(np.asarray(axis_ranges[sidx]), i))
                else:
                    raise ValueError(incompatible_axis_ranges_error_message(i))

        return axis_ranges, extent

    @abc.abstractmethod
    def plot_start_image(self, ax):
        """
        Abstract method for plotting first slice of array.

        Must exist here but be defined in subclass.
        """

    @abc.abstractmethod
    def update_plot(self, val, artist, slider):
        """
        Abstract method for updating the plot.

        Must exist here but be defined in subclass.
        """
        ind = int(val)
        ax_ind = self.slider_axes[slider.slider_ind]
        # Update slider label to reflect real world values in axis_ranges.
        label = self.axis_ranges[ax_ind](ind)
        if isinstance(label, u.Quantity):
            slider.valtext.set_text(label.to_string(precision=5,
                                                    format='latex',
                                                    subfmt='inline'))
        elif isinstance(label, str):
            slider.valtext.set_text(label)
        else:
            slider.valtext.set_text(f"{label:10.2f}")


def edges_to_centers_nd(axis_range, edges_axis):
    """
    Converts ND array of pixel edges to pixel centers along one axis.

    Parameters
    ----------
    axis_range: `numpy.ndarray`
        Array of pixel edges.
    edges_axis: `int`
        Index of axis along which centers are to be calculated.
    """
    upper_edge_indices = [slice(None)] * axis_range.ndim
    upper_edge_indices[edges_axis] = slice(1, axis_range.shape[edges_axis])
    upper_edges = axis_range[tuple(upper_edge_indices)]

    lower_edge_indices = [slice(None)] * axis_range.ndim
    lower_edge_indices[edges_axis] = slice(0, -1)
    lower_edges = axis_range[tuple(lower_edge_indices)]

    return (upper_edges - lower_edges) / 2 + lower_edges
