"""
This module provides plotting support in iPython.
"""
from functools import wraps
from contextlib import contextmanager

import matplotlib.pyplot as plt

from astropy.coordinates import BaseCoordinateFrame, SkyCoord

from sunpy.coordinates import RotatedSunFrame
from .wcsaxes_compat import is_wcsaxes

__all__ = ['peek_show', "axis_labels_from_ctype", 'diffrot_axes']


def peek_show(func):
    """
    A decorator to place on ``peek()`` methods to show the figure.

    The ``peek()`` method should return the figure then this method will
    attempt to show it in the correct way. This decorator will not return the
    figure to the user.
    """
    @wraps(func)
    def show_figure(*args, **kwargs):
        _ = func(*args, **kwargs)
        plt.show()

    return show_figure


def axis_labels_from_ctype(ctype, unit):
    """
    Returns axis labels for the given coordinate type and unit.

    Parameters
    ----------
    ctype: `str`
        Coordinate type.
    unit: `str`
        Required unit.

    Returns
    -------
    `str`
        "Axis Label [Unit]"
    """
    ctype_short = ctype[:4]

    labels = {'HGLN': f'Heliographic Longitude [{unit}]',
              'CRLN': f'Carrington Longitude [{unit}]',
              'HPLN': f'Helioprojective Longitude (Solar-X) [{unit}]',
              'SOLX': f'Heliocentric X [{unit}]',

              'HGLT': f'Latitude [{unit}]',
              'CRLT': f'Latitude [{unit}]',
              'HPLT': f'Helioprojective Latitude (Solar-Y) [{unit}]',
              'SOLY': f'Heliocentric Y [{unit}]'}

    return labels.get(ctype_short, f"{ctype} [{unit}]")


@contextmanager
def diffrot_axes(axes=None):
    """
    Context manager that enables a WCSAxes instance to automatically account for
    solar differential rotation.

    Coordinates and coordinate frames supplied to WCSAxes methods (e.g.,
    :meth:`~astropy.visualization.wcsaxes.WCSAxes.plot_coord`,
    :meth:`~astropy.visualization.wcsaxes.WCSAxes.get_coords_overlay`,
    :meth:`~astropy.visualization.wcsaxes.WCSAxes.get_transform`) account for the
    differential rotation between their ``obstime`` and the ``obstime`` of the
    WCSAxes's coordinate frame.

    Parameters
    ----------
    axes : `~astropy.visualization.wcsaxes.WCSAxes`
         The axes to modify.  If ``None``, the current `matplotlib.pyplot` axes are
         used.

    Notes
    -----
    Since this context manager modifies WCSAxes methods, it may not function as
    intended if the WCSAxes implementation changes.

    This context manager uses `~sunpy.coordinates.metaframes.RotatedSunFrame`
    internally.
    """
    try:
        if axes is None:
            axes = plt.gca()

        if not is_wcsaxes(axes):
            raise TypeError("The axes need to be an instance of WCSAxes. You may have neglected "
                            "to use `projection=` when creating the axes.")

        # Monkey patch the instance's private method used by get_transform()
        old_get_transform = axes._get_transform_no_transdata

        def new_get_transform(self, frame):
            transform = old_get_transform(frame)
            if getattr(transform, '_a', None) is axes._transform_pixel2world:
                step = transform._b
                if step.output_system.obstime is None:
                    raise ValueError(f"The {step.output_system.__class__.__name__} instance "
                                     f"must have a defined `obstime`.")
                new_input_system = RotatedSunFrame(base=step.input_system,
                                                   rotated_time=step.output_system.obstime)
                transform._b = step.__class__(new_input_system, step.output_system)
            return transform

        axes._get_transform_no_transdata = new_get_transform.__get__(axes, type(axes))

        # Monkey patch the instance's plot_coord()
        # This is necessary because plot_coord() first transforms the coordinate prior to calling
        #   get_transform()
        old_plot_coord = axes.plot_coord

        def new_plot_coord(self, *args, **kwargs):
            if isinstance(args[0], (SkyCoord, BaseCoordinateFrame)):
                frame = args[0].frame if hasattr(args[0], 'frame') else args[0]
                if frame.obstime is None:
                    raise ValueError(f"The {frame.__class__.__name__} instance "
                                     f"must have a defined `obstime`.")
                rsf_frame = RotatedSunFrame(base=axes._transform_pixel2world.frame_out,
                                            rotated_time=frame.obstime)
                args = (frame.transform_to(rsf_frame).as_base(),) + args[1:]
            return old_plot_coord(*args, **kwargs)

        axes.plot_coord = new_plot_coord.__get__(axes, type(axes))

        yield

    finally:
        # Reverse the monkey patches
        axes._get_transform_no_transdata = old_get_transform
        axes.plot_coord = old_plot_coord
