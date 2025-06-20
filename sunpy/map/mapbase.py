"""
Map is a generic Map class from which all other Map classes inherit from.
"""
import re
import copy
import html
import inspect
import numbers
import textwrap
import warnings
import itertools
import webbrowser
from typing import Literal
from tempfile import NamedTemporaryFile
from collections import namedtuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import reproject
from matplotlib.backend_bases import FigureCanvasBase
from matplotlib.figure import Figure

try:
    from dask.array import Array as DaskArray
    DASK_INSTALLED = True
except ImportError:
    DASK_INSTALLED = False

import astropy.units as u
import astropy.wcs
from astropy.coordinates import BaseCoordinateFrame, Longitude, SkyCoord, UnitSphericalRepresentation
from astropy.nddata import NDData
from astropy.utils.metadata import MetaData
from astropy.visualization import HistEqStretch, ImageNormalize
from astropy.visualization.wcsaxes import Quadrangle, WCSAxes

# The next two are not used but are called to register functions with external modules
import sunpy.coordinates
import sunpy.visualization.colormaps
from sunpy import config, log
from sunpy.coordinates import HeliographicCarrington, get_earth, sun
from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.image.resample import resample as sunpy_image_resample
from sunpy.image.resample import reshape_image_to_4d_superpixel
from sunpy.image.transform import _get_transform_method, _rotation_function_names, affine_transform
from sunpy.io._file_tools import write_file
from sunpy.io._fits import extract_waveunit, header_to_fits
from sunpy.map.maputils import _clip_interval, _handle_norm
from sunpy.sun import constants
from sunpy.time import is_time, parse_time
from sunpy.util import MetaDict, expand_list, extent_in_other_wcs, grid_perimeter
from sunpy.util.decorators import (
    add_common_docstring,
    cached_property_based_on,
    check_arithmetic_compatibility,
    deprecated,
)
from sunpy.util.exceptions import SunpyUserWarning, warn_deprecated, warn_metadata, warn_user
from sunpy.util.functools import seconddispatch
from sunpy.util.util import _figure_to_base64, fix_duplicate_notes
from sunpy.visualization import axis_labels_from_ctype, peek_show, wcsaxes_compat
from sunpy.visualization.colormaps import cm as sunpy_cm
from sunpy.visualization.visualization import _PrecomputedPixelCornersTransform

TIME_FORMAT = config.get("general", "time_format")
PixelPair = namedtuple('PixelPair', 'x y')
SpatialPair = namedtuple('SpatialPair', 'axis1 axis2')
_NUMPY_COPY_IF_NEEDED = False if np.__version__.startswith("1.") else None
_META_FIX_URL = 'https://docs.sunpy.org/en/stable/how_to/fix_map_metadata.html'

# Manually specify the ``.meta`` docstring. This is assigned to the .meta
# class attribute in GenericMap.__init__()
_meta_doc = """
The map metadata.

This is used to interpret the map data. It may
have been modified from the original metadata by sunpy. See the
`~sunpy.util.MetaDict.added_items`, `~sunpy.util.MetaDict.removed_items`
and `~sunpy.util.MetaDict.modified_items` properties of MetaDict
to query how the metadata has been modified.
"""

# The notes live here so we can reuse it in the source maps
_notes_doc = """

Notes
-----

A number of the properties of this class are returned as two-value named
tuples that can either be indexed by position ([0] or [1]) or be accessed
by the names (.x and .y) or (.axis1 and .axis2). Things that refer to pixel
axes use the ``.x``, ``.y`` convention, where x and y refer to the FITS
axes (x for columns y for rows). Spatial axes use ``.axis1`` and ``.axis2``
which correspond to the first and second axes in the header. ``axis1``
corresponds to the coordinate axis for ``x`` and ``axis2`` corresponds to
``y``.

This class assumes that the metadata adheres to the FITS 4 standard.
Where the CROTA2 metadata is provided (without PC_ij) it assumes a conversion
to the standard PC_ij described in section 6.1 of :cite:t:`calabretta_representations_2002`.

.. warning::
    If a header has CD_ij values but no PC_ij values, CDELT values are required
    for this class to construct the WCS.
    If a file with more than two dimensions is feed into the class,
    only the first two dimensions (NAXIS1, NAXIS2) will be loaded and the
    rest will be discarded.
"""

__all__ = ['GenericMap', 'MapMetaValidationError', 'PixelPair']


class MapMetaValidationError(AttributeError):
    pass


class GenericMap(NDData):
    """
    A Generic spatially-aware 2D data array

    Parameters
    ----------
    data : `numpy.ndarray`, list
        A 2d list or ndarray containing the map data.
    header : dict
        A dictionary of the original image header tags.
    plot_settings : dict, optional
        Plot settings.

    Other Parameters
    ----------------
    **kwargs :
        Additional keyword arguments are passed to `~astropy.nddata.NDData`
        init.


    Methods and their known behavior with dask arrays
    -------------------------------------------------

    +-------------------+------------------------------------+
    | Method            | Preserve laziness with Dask Arrays |
    +===================+====================================+
    | `reproject_to`    | No                                 |
    +-------------------+------------------------------------+
    | `resample`        | No                                 |
    +-------------------+------------------------------------+
    | `rotate`          | No                                 |
    +-------------------+------------------------------------+
    | `max`             | Yes                                |
    +-------------------+------------------------------------+
    | `mean`            | Yes                                |
    +-------------------+------------------------------------+
    | `min`             | Yes                                |
    +-------------------+------------------------------------+
    | `std`             | Yes                                |
    +-------------------+------------------------------------+
    | `superpixel`      | Yes                                |
    +-------------------+------------------------------------+
    | `submap`          | Yes                                |
    +-------------------+------------------------------------+


    Examples
    --------
    >>> import sunpy.map
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
    >>> aia   # doctest: +REMOTE_DATA
    <sunpy.map.sources.sdo.AIAMap object at ...>
    SunPy Map
    ---------
    Observatory:                 SDO
    Instrument:          AIA 3
    Detector:            AIA
    Measurement:                 171.0 Angstrom
    Wavelength:          171.0 Angstrom
    Observation Date:    2011-06-07 06:33:02
    Exposure Time:               0.234256 s
    Dimension:           [1024. 1024.] pix
    Coordinate System:   helioprojective
    Scale:                       [2.402792 2.402792] arcsec / pix
    Reference Pixel:     [511.5 511.5] pix
    Reference Coord:     [3.22309951 1.38578135] arcsec
    array([[ -95.92475  ,    7.076416 ,   -1.9656711, ..., -127.96519  ,
            -127.96519  , -127.96519  ],
           [ -96.97533  ,   -5.1167884,    0.       , ...,  -98.924576 ,
            -104.04137  , -127.919716 ],
           [ -93.99607  ,    1.0189276,   -4.0757103, ...,   -5.094638 ,
             -37.95505  , -127.87541  ],
           ...,
           [-128.01454  , -128.01454  , -128.01454  , ..., -128.01454  ,
            -128.01454  , -128.01454  ],
           [-127.899666 , -127.899666 , -127.899666 , ..., -127.899666 ,
            -127.899666 , -127.899666 ],
           [-128.03072  , -128.03072  , -128.03072  , ..., -128.03072  ,
            -128.03072  , -128.03072  ]], shape=(1024, 1024), dtype=float32)

    >>> aia.spatial_units   # doctest: +REMOTE_DATA
    SpatialPair(axis1=Unit("arcsec"), axis2=Unit("arcsec"))
    >>> aia.peek()   # doctest: +SKIP
    """
    _registry = dict()
    # This overrides the default doc for the meta attribute
    meta = MetaData(doc=_meta_doc, copy=False)
    # Enabling the GenericMap reflected operators is a bit subtle. The GenericMap
    # reflected operator will be used only if the Quantity non-reflected operator
    # returns NotImplemented. The Quantity operator strips the unit from the
    # Quantity and tries to combine the value with the GenericMap using NumPy's
    # __array_ufunc__(). If NumPy believes that it can proceed, this will result
    # in an error. We explicitly set __array_ufunc__ = None so that the NumPy
    # call, and consequently the Quantity operator, will return NotImplemented.
    __array_ufunc__ = None

    def __init_subclass__(cls, **kwargs):
        """
        An __init_subclass__ hook initializes all of the subclasses of a given class.
        So for each subclass, it will call this block of code on import.
        This replicates some metaclass magic without the need to be aware of metaclasses.
        Here we use this to register each subclass in a dict that has the
        ``is_datasource_for`` attribute.
        This is then passed into the Map Factory so we can register them.
        """
        super().__init_subclass__(**kwargs)
        if cls.__doc__ is None:
            # Set an empty string, to prevent an error adding None to str in the next line
            cls.__doc__ = ''
        cls.__doc__ = fix_duplicate_notes(_notes_doc, cls.__doc__)

        if hasattr(cls, 'is_datasource_for'):
            # NOTE: This conditional is due to overlapping map sources in sunpy and pfsspy that
            # lead to a MultipleMatchError if sunpy.map and pfsspy.map are imported.
            # See https://github.com/sunpy/sunpy/issues/7294 for more information.
            # This also applies to older versions of sunkit-magex with ADAPTMap.
            if f'{cls.__module__}.{cls.__name__}' not in  ["pfsspy.map.GongSynopticMap", "sunkit_magex.pfss.map.ADAPTMap"]:
                cls._registry[cls] = cls.is_datasource_for

    def __init__(self, data, header, plot_settings=None, **kwargs):
        # If the data has more than two dimensions, the first dimensions
        # (NAXIS1, NAXIS2) are used and the rest are discarded.
        ndim = data.ndim
        if ndim > 2:
            # We create a slice that removes all but the 'last' two
            # dimensions. (Note dimensions in ndarray are in reverse order)

            new_2d_slice = [0]*(ndim-2)
            new_2d_slice.extend([slice(None), slice(None)])
            data = data[tuple(new_2d_slice)]
            # Warn the user that the data has been truncated
            warn_user("This file contains more than 2 dimensions. "
                      "Data will be truncated to the first two dimensions.")

        params = list(inspect.signature(NDData).parameters)
        nddata_kwargs = {x: kwargs.pop(x) for x in params & kwargs.keys()}
        super().__init__(data, meta=MetaDict(header), **nddata_kwargs)

        # Setup some attributes
        self._nickname = None
        # These are placeholders for default attributes, which are only set
        # once if their data isn't present in the map metadata.
        self._default_time = None
        self._default_dsun = None
        self._default_carrington_longitude = None
        self._default_heliographic_latitude = None
        self._default_heliographic_longitude = None

        # Validate header
        # TODO: This should be a function of the header, not of the map
        self._validate_meta()
        self.plot_settings = {'cmap': 'gray',
                            'interpolation': 'nearest',
                            'origin': 'lower'
                            }
        if self.dtype != np.uint8:
            # Put import here to reduce sunpy.map import time
            from matplotlib import colors
            self.plot_settings['norm'] = colors.Normalize()
        if plot_settings:
            self.plot_settings.update(plot_settings)

        # Try and set the colormap. This is not always possible if this method
        # is run before map sources fix some of their metadata, so
        # just ignore any exceptions raised.
        try:
            cmap = self._get_cmap_name()
            if cmap in sunpy_cm.cmlist:
                self.plot_settings['cmap'] = cmap
        except Exception:
            pass

    def __getitem__(self, key):
        """ This should allow indexing by physical coordinate """
        raise NotImplementedError(
            "The ability to index Map by physical"
            " coordinate is not yet implemented.")

    def _text_summary(self):
        dt = self.exposure_time
        wave = self.wavelength
        measurement = self.measurement

        dt = 'Unknown' if dt is None else dt
        wave = 'Unknown' if wave is None else wave
        measurement = 'Unknown' if measurement is None else measurement

        return textwrap.dedent("""\
                   SunPy Map
                   ---------
                   Observatory:\t\t {obs}
                   Instrument:\t\t {inst}
                   Detector:\t\t {det}
                   Measurement:\t\t {meas}
                   Wavelength:\t\t {wave}
                   Observation Date:\t {date}
                   Exposure Time:\t\t {dt}
                   Dimension:\t\t {dim}
                   Coordinate System:\t {coord}
                   Scale:\t\t\t {scale}
                   Reference Pixel:\t {refpix}
                   Reference Coord:\t {refcoord}\
                   """).format(obs=self.observatory, inst=self.instrument, det=self.detector,
                               meas=measurement, wave=wave,
                               date=self.date.strftime(TIME_FORMAT),
                               dt=dt,
                               dim=u.Quantity(self.dimensions),
                               scale=u.Quantity(self.scale),
                               coord=self._coordinate_frame_name,
                               refpix=u.Quantity(self.reference_pixel),
                               refcoord=u.Quantity((self._reference_longitude,
                                                    self._reference_latitude)),
                               tmf=TIME_FORMAT)

    def __str__(self):
        return f"{self._text_summary()}\n{self.data.__repr__()}"

    def __repr__(self):
        return f"{object.__repr__(self)}\n{self}"

    def _repr_html_(self, compute_dask=False):
        """
        Produce an HTML summary with plots for use in Jupyter notebooks.
        """
        # Convert the text repr to an HTML table
        partial_html = self._text_summary()[20:].replace('\n', '</td></tr><tr><th>')\
                                                .replace(':\t', '</th><td>')
        text_to_table = textwrap.dedent(f"""\
            <table style='text-align:left'>
                <tr><th>{partial_html}</td></tr>
            </table>""").replace('\n', '')

        # Handle bad values (infinite and NaN) in the data array
        finite_data = self.data[np.isfinite(self.data)]
        count_nan = np.isnan(self.data).sum()
        count_inf = np.isinf(self.data).sum()

        if DASK_INSTALLED and isinstance(finite_data, DaskArray):
            # This will fetch the entire data array into memory and only happens for the quicklook method
            if compute_dask:
                finite_data = finite_data.compute()
            else:
                dask_html = self.data._repr_html_()
                return textwrap.dedent(f"""\
                    <pre>{html.escape(object.__repr__(self))}</pre>
                    <table>
                        <tr>
                            <td>{text_to_table}</td>
                            <td>
                                {dask_html}
                            </td>
                        </tr>
                        <tr>
                        </tr>
                    </table>""")

        # Assemble an informational string with the counts of bad pixels
        bad_pixel_text = ""
        if count_nan + count_inf > 0:
            bad_pixel_text = "Bad pixels are shown in red: "
            text_list = []
            if count_nan > 0:
                text_list.append(f"{count_nan} NaN")
            if count_inf > 0:
                text_list.append(f"{count_inf} infinite")
            bad_pixel_text += ", ".join(text_list)

        # Use a grayscale colormap with histogram equalization (and red for bad values)
        # Make a copy of the colormap to avoid modifying the matplotlib instance when
        # doing set_bad() (copy not needed when min mpl is 3.5, as already a copy)
        cmap = copy.copy(matplotlib.colormaps['gray'])
        cmap.set_bad(color='red')
        norm = ImageNormalize(stretch=HistEqStretch(finite_data))

        # Plot the image in pixel space
        fig = Figure(figsize=(5.2, 4.8))
        # Figure instances in matplotlib<3.1 do not create a canvas by default
        if fig.canvas is None:
            FigureCanvasBase(fig)
        ax = fig.subplots()
        ax.imshow(self.data, origin='lower', interpolation='nearest', cmap=cmap, norm=norm)
        ax.set_xlabel('X pixel')
        ax.set_ylabel('Y pixel')
        ax.set_title('In pixel space')
        pixel_src = _figure_to_base64(fig)
        bounds = ax.get_position().bounds  # save these axes bounds for later use

        # Plot the image using WCS information, with the same axes bounds as above
        fig = Figure(figsize=(5.2, 4.8))
        # Figure instances in matplotlib<3.1 do not create a canvas by default
        if fig.canvas is None:
            FigureCanvasBase(fig)
        # Create the WCSAxes manually because we need to avoid using pyplot
        ax = WCSAxes(fig, bounds, aspect='equal', wcs=self.wcs)
        fig.add_axes(ax)
        self.plot(axes=ax, cmap=cmap, norm=norm)
        ax.set_title('In coordinate space using WCS information')
        wcs_src = _figure_to_base64(fig)

        # Plot the histogram of pixel values
        fig = Figure(figsize=(4.8, 2.4), constrained_layout=True)
        # Figure instances in matplotlib<3.1 do not create a canvas by default
        if fig.canvas is None:
            FigureCanvasBase(fig)
        ax = fig.subplots()
        values, bins, patches = ax.hist(finite_data.ravel(), bins=100)
        norm_centers = norm(0.5 * (bins[:-1] + bins[1:])).data
        for c, p in zip(norm_centers, patches):
            plt.setp(p, "facecolor", cmap(c))
        ax.plot(np.array([bins[:-1], bins[1:]]).T.ravel(),
                np.array([values, values]).T.ravel())
        ax.set_facecolor('white')
        ax.semilogy()
        # Explicitly set the power limits for the X axis formatter to avoid text overlaps
        ax.xaxis.get_major_formatter().set_powerlimits((-3, 4))
        ax.set_xlabel(f"Pixel value{' (' + str(self.unit) + ')' if self.unit else ''} in linear bins")
        ax.set_ylabel('# of pixels')
        ax.set_title('Distribution of pixel values [click for cumulative]')
        hist_src = _figure_to_base64(fig)

        # Plot the CDF of the pixel values using a symmetric-log horizontal scale
        fig = Figure(figsize=(4.8, 2.4), constrained_layout=True)
        # TODO: Figure instances in matplotlib<3.1 do not create a canvas by default
        if fig.canvas is None:
            FigureCanvasBase(fig)
        ax = fig.subplots()
        n_bins = 256
        bins = norm.inverse(np.arange(n_bins + 1) / n_bins)
        values, _, patches = ax.hist(finite_data.ravel(), bins=bins, cumulative=True)
        for i, p in enumerate(patches):
            plt.setp(p, "facecolor", cmap((i + 0.5) / n_bins))
        ax.plot(np.array([bins[:-1], bins[1:]]).T.ravel(),
                np.array([values, values]).T.ravel())
        ax.set_facecolor('white')
        ax.set_xscale('symlog')
        ax.set_yscale('log')
        ax.set_xlabel(f"Pixel value{' (' + str(self.unit) + ')' if self.unit else ''} in equalized bins")
        ax.set_ylabel('Cumulative # of pixels')
        ax.set_title('Cumulative distribution of pixel values')
        cdf_src = _figure_to_base64(fig)

        return textwrap.dedent(f"""\
            <pre>{html.escape(object.__repr__(self))}</pre>
            <table>
                <tr>
                    <td>{text_to_table}</td>
                    <td rowspan=3>
                        <div align=center>
                            Image colormap uses histogram equalization<br>
                            Click on the image to toggle between units
                        </div>
                        <img src='data:image/png;base64,{wcs_src}'
                             src2='data:image/png;base64,{pixel_src}'
                             onClick='var temp = this.src;
                                      this.src = this.getAttribute("src2");
                                      this.setAttribute("src2", temp)'
                        />
                        <div align=center>
                            {bad_pixel_text}
                        </div>
                    </td>
                </tr>
                <tr>
                </tr>
                <tr>
                    <td><img src='data:image/png;base64,{hist_src}'
                             src2='data:image/png;base64,{cdf_src}'
                             onClick='var temp = this.src;
                                      this.src = this.getAttribute("src2");
                                      this.setAttribute("src2", temp)'
                        />
                    </td>
                </tr>
            </table>""")

    def quicklook(self):
        """
        Display a quicklook summary of the Map instance using the default web browser.

        Notes
        -----
        The image colormap uses
        `histogram equalization <https://en.wikipedia.org/wiki/Histogram_equalization>`__.

        Clicking on the image to switch between pixel space and coordinate space requires
        Javascript support to be enabled in the web browser.

        Examples
        --------
        >>> from sunpy.map import Map
        >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
        >>> smap = Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
        >>> smap.quicklook()  # doctest: +SKIP

        (which will open the following content in the default web browser)

        .. generate:: html
            :html_border:

            from sunpy.map import Map
            import sunpy.data.sample
            smap = Map(sunpy.data.sample.AIA_171_IMAGE)
            print(smap._repr_html_())

        """
        with NamedTemporaryFile('w', delete=False, prefix='sunpy.map.', suffix='.html') as f:
            url = 'file://' + f.name
            f.write(textwrap.dedent(f"""\
                <html>
                    <title>Quicklook summary for {html.escape(object.__repr__(self))}</title>
                    <body>{self._repr_html_(compute_dask=True)}</body>
                </html>"""))
        webbrowser.open_new_tab(url)

    @classmethod
    def _new_instance(cls, data, meta, plot_settings=None, **kwargs):
        """
        Instantiate a new instance of this class using given data.
        This is a shortcut for ``type(self)(data, meta, plot_settings)``.
        """
        new_map = cls(data, meta, **kwargs)
        # plot_settings are set explicitly here as some map sources
        # explicitly set some of the plot_settings in the constructor
        # and we want to preserve the plot_settings of the previous
        # instance.
        if plot_settings is not None:
            new_map.plot_settings.update(plot_settings)
        return new_map

    def _get_lon_lat(self, frame):
        """
        Given a coordinate frame, extract the lon and lat by casting to
        SphericalRepresentation first.
        """
        r = frame.represent_as(UnitSphericalRepresentation)
        return r.lon.to(self.spatial_units[0]), r.lat.to(self.spatial_units[1])

    @property
    def quantity(self):
        """Unitful representation of the map data."""
        return u.Quantity(self.data, self.unit, copy=_NUMPY_COPY_IF_NEEDED)

    def _new_instance_from_op(self, new_data):
        """
        Helper function for creating new map instances after arithmetic
        operations.
        """
        new_meta = copy.deepcopy(self.meta)
        new_meta['bunit'] = new_data.unit.to_string()
        return self._new_instance(new_data.value, new_meta, plot_settings=self.plot_settings)

    def __neg__(self):
        return self._new_instance(-self.data, self.meta, plot_settings=self.plot_settings)

    @check_arithmetic_compatibility
    def __pow__(self, value):
        new_data = self.quantity ** value
        return self._new_instance_from_op(new_data)

    @check_arithmetic_compatibility
    def __add__(self, value):
        new_data = self.quantity + value
        return self._new_instance_from_op(new_data)

    def __radd__(self, value):
        return self.__add__(value)

    def __sub__(self, value):
        return self.__add__(-value)

    def __rsub__(self, value):
        return self.__neg__().__add__(value)

    @check_arithmetic_compatibility
    def __mul__(self, value):
        new_data = self.quantity * value
        return self._new_instance_from_op(new_data)

    def __rmul__(self, value):
        return self.__mul__(value)

    @check_arithmetic_compatibility
    def __truediv__(self, value):
        return self.__mul__(1/value)

    @check_arithmetic_compatibility
    def __rtruediv__(self, value):
        new_data = value / self.quantity
        return self._new_instance_from_op(new_data)

    @property
    def _meta_hash(self):
        return self.meta.item_hash()

    def _set_symmetric_vmin_vmax(self):
        """
        Set symmetric vmin and vmax about zero
        """
        threshold = np.nanmax(abs(self.data))
        self.plot_settings['norm'].vmin = -threshold
        self.plot_settings['norm'].vmax = threshold

    @property
    @cached_property_based_on('_meta_hash')
    def wcs(self):
        """
        The `~astropy.wcs.WCS` property of the map.

        Notes
        -----
        ``dateobs`` is always populated with the "canonical" observation time as
        provided by the `.date` property. This will commonly be the DATE-OBS key if it
        is in the metadata, but see that property for the logic otherwise.

        ``dateavg`` is always populated with the reference date of the coordinate system
        as provided by the `.reference_date` property. This will commonly be the
        DATE-AVG key if it is in the metadata, but see that property for the logic
        otherwise.

        ``datebeg`` is conditonally populated with the start of the observation period
        as provided by the `.date_start` property, which normally returns a value only
        if the DATE-BEG key is in the metadata.

        ``dateend`` is conditonally populated with the end of the observation period as
        provided by the `.date_end` property, which normally returns a value only if the
        DATE-END key is in the metadata.
        """
        w2 = astropy.wcs.WCS(naxis=2)

        # Add one to go from zero-based to one-based indexing
        w2.wcs.crpix = u.Quantity(self.reference_pixel) + 1 * u.pix
        # Make these a quantity array to prevent the numpy setting element of
        # array with sequence error.
        # Explicitly call ``.to()`` to check that scale is in the correct units
        w2.wcs.cdelt = u.Quantity([self.scale[0].to(self.spatial_units[0] / u.pix),
                                   self.scale[1].to(self.spatial_units[1] / u.pix)])
        w2.wcs.crval = u.Quantity([self._reference_longitude, self._reference_latitude])
        w2.wcs.ctype = self.coordinate_system
        w2.wcs.pc = self.rotation_matrix
        w2.wcs.set_pv(self._pv_values)
        # FITS standard doesn't allow both PC_ij *and* CROTA keywords
        w2.wcs.crota = (0, 0)
        w2.wcs.cunit = self.spatial_units
        w2.wcs.aux.rsun_ref = self.rsun_meters.to_value(u.m)

        w2.wcs.dateobs = self.date.utc.isot
        w2.wcs.dateavg = self.reference_date.utc.isot
        if self.date_start is not None:
            w2.wcs.datebeg = self.date_start.utc.isot
        if self.date_end is not None:
            w2.wcs.dateend = self.date_end.utc.isot

        # Set observer coordinate information except when we know it is not appropriate (e.g., HGS)
        sunpy_frame = sunpy.coordinates.wcs_utils._sunpy_frame_class_from_ctypes(w2.wcs.ctype)
        if sunpy_frame is None or hasattr(sunpy_frame, 'observer'):
            # Clear all the aux information that was set earlier. This is to avoid
            # issues with maps that store multiple observer coordinate keywords.
            # Note that we have to create a new WCS as it's not possible to modify
            # wcs.wcs.aux in place.
            header = w2.to_header()
            for kw in ['crln_obs', 'dsun_obs', 'hgln_obs', 'hglt_obs']:
                header.pop(kw, None)
            w2 = astropy.wcs.WCS(header)

            # Get observer coord, and set the aux information
            obs_coord = self.observer_coordinate
            sunpy.coordinates.wcs_utils._set_wcs_aux_obs_coord(w2, obs_coord)

        # Set the shape of the data array
        w2.array_shape = self.data.shape

        # Validate the WCS here.
        w2.wcs.set()
        return w2

    @property
    def coordinate_frame(self):
        """
        An `astropy.coordinates.BaseCoordinateFrame` instance created from the coordinate
        information for this Map, or None if the frame cannot be determined.

        Notes
        -----
        The ``obstime`` for the coordinate frame uses the `.reference_date` property,
        which may be different from the `.date` property.
        """
        try:
            return astropy.wcs.utils.wcs_to_celestial_frame(self.wcs)
        except ValueError as e:
            warn_user(f'Could not determine coordinate frame from map metadata.\n{e}')
            return None

    @property
    def _coordinate_frame_name(self):
        if self.coordinate_frame is None:
            return 'Unknown'
        return self.coordinate_frame.name

    def _as_mpl_axes(self):
        """
        Compatibility hook for Matplotlib and WCSAxes.
        This functionality requires the WCSAxes package to work. The reason
        we include this here is that it allows users to use WCSAxes without
        having to explicitly import WCSAxes
        With this method, one can do::

            import matplotlib.pyplot as plt
            import sunpy.map
            amap = sunpy.map.Map('filename.fits')
            fig = plt.figure()
            ax = plt.subplot(projection=amap)
            ...

        and this will generate a plot with the correct WCS coordinates on the
        axes. See <https://docs.astropy.org/en/stable/visualization/wcsaxes/index.html> for more information.
        """
        # This code is reused from Astropy
        return WCSAxes, {'wcs': self.wcs}

    # Some numpy extraction
    @property
    def dimensions(self):
        """
        The dimensions of the array (x axis first, y axis second).
        """
        return PixelPair(*u.Quantity(np.flipud(self.data.shape), 'pixel'))

    @property
    def dtype(self):
        """
        The `numpy.dtype` of the array of the map.
        """
        return self.data.dtype

    @property
    def ndim(self):
        """
        The value of `numpy.ndarray.ndim` of the data array of the map.
        """
        return self.data.ndim

    def std(self, *args, **kwargs):
        """
        Calculate the standard deviation of the data array, ignoring NaNs.

        This method **does** preserve dask arrays.
        """
        return np.nanstd(self.data, *args, **kwargs)

    def mean(self, *args, **kwargs):
        """
        Calculate the mean of the data array, ignoring NaNs.

        This method **does** preserve dask arrays.
        """
        return np.nanmean(self.data, *args, **kwargs)

    def min(self, *args, **kwargs):
        """
        Calculate the minimum value of the data array, ignoring NaNs.

        This method **does** preserve dask arrays.
        """
        return np.nanmin(self.data, *args, **kwargs)

    def max(self, *args, **kwargs):
        """
        Calculate the maximum value of the data array, ignoring NaNs.

        This method **does** preserve dask arrays.
        """
        return np.nanmax(self.data, *args, **kwargs)

    @staticmethod
    def _parse_fits_unit(unit_str):
        replacements = {'gauss': 'G',
                        'counts / pixel': 'ct/pix',}
        if unit_str.lower() in replacements:
            unit_str = replacements[unit_str.lower()]
        unit = u.Unit(unit_str, parse_strict='silent')
        for base in unit.bases:
            # NOTE: Special case DN here as it is not part of the FITS standard, but
            # is widely used and is also a recognized astropy unit
            if base is u.DN:
                continue
            try:
                if isinstance(base, u.UnrecognizedUnit):
                    raise ValueError

                # Also rejects a unit that is not in the FITS standard but is equivalent to one (e.g., Mx)
                if u.Unit(base.to_string(format='fits')) is not base:  # to_string() can raise ValueError
                    raise ValueError
            except ValueError:
                warn_metadata(f'Could not parse unit string "{unit_str}" as a valid FITS unit.\n'
                              f'See {_META_FIX_URL} for how to fix metadata before loading it '
                               'with sunpy.map.Map.\n'
                               'See https://fits.gsfc.nasa.gov/fits_standard.html for '
                               'the FITS unit standards.')
                return None
        return unit

    @property
    def unit(self):
        """
        Unit of the map data.

        This is taken from the 'BUNIT' FITS keyword. If no 'BUNIT' entry is
        present in the metadata then this returns `None`. If the 'BUNIT' value
        cannot be parsed into a unit a warning is raised, and `None` returned.
        """
        unit_str = self.meta.get('bunit', None)
        if unit_str is None:
            return
        return self._parse_fits_unit(unit_str)

# #### Keyword attribute and other attribute definitions #### #

    def _base_name(self):
        """Abstract the shared bit between name and latex_name"""
        if self.measurement is None:
            format_str = "{nickname} {date}"
        else:
            format_str = "{nickname} {{measurement}} {date}"
        return format_str.format(nickname=self.nickname,
                                 date=parse_time(self.date).strftime(TIME_FORMAT))

    @property
    def name(self):
        """Human-readable description of the Map."""
        return self._base_name().format(measurement=self.measurement)

    @property
    def latex_name(self):
        """LaTeX formatted description of the Map."""
        if isinstance(self.measurement, u.Quantity):
            return self._base_name().format(measurement=self.measurement._repr_latex_())
        else:
            return self.name

    @property
    def nickname(self):
        """An abbreviated human-readable description of the map-type; part of
        the Helioviewer data model."""
        return self._nickname if self._nickname else self.detector

    @nickname.setter
    def nickname(self, n):
        self._nickname = n

    def _get_date(self, key):
        time = self.meta.get(key, None)
        if not time:
            return

        # Get the time scale
        if 'TAI' in time:
            # SDO specifies the 'TAI' scale in their time string, which is parsed
            # by parse_time(). If a different timescale is also present, warn the
            # user that it will be ignored.
            timesys = 'TAI'
            timesys_meta = self.meta.get('timesys', '').upper()
            if timesys_meta not in ('', 'TAI'):
                warn_metadata('Found "TAI" in time string, ignoring TIMESYS keyword '
                              f'which is set to "{timesys_meta}".')
        else:
            timesys = self._timesys

        return parse_time(time, scale=timesys.lower())

    @property
    def _timesys(self):
        """
        Time system.
        """
        # UTC is the FITS standard default
        return self.meta.get('timesys', 'UTC')

    @property
    def date_start(self):
        """
        Time of the beginning of the image acquisition.

        Taken from the DATE-BEG FITS keyword.
        """
        return self._get_date('date-beg')

    @property
    def date_end(self):
        """
        Time of the end of the image acquisition.

        Taken from the DATE-END FITS keyword.
        """
        return self._get_date('date-end')

    @property
    def date_average(self):
        """
        Average time of the image acquisition.

        Taken from the DATE-AVG FITS keyword if present, otherwise halfway
        between `date_start` and `date_end` if both pieces of metadata are
        present.
        """
        avg = self._get_date('date-avg')
        if avg is None:
            start, end = self.date_start, self.date_end
            if start is not None and end is not None:
                avg = start + (end - start) / 2

        return avg

    @property
    def _date_obs(self):
        # Get observation date from date-obs, falling back to date_obs
        if is_time(self.meta.get("date-obs", None)):
            return self._get_date('date-obs')
        elif is_time(self.meta.get('date_obs', None)):
            return self._get_date('date_obs')

    @property
    def reference_date(self):
        """
        The reference date for the coordinate system.

        This date is used to define the ``obstime`` of the coordinate frame and often
        the ``obstime`` of the observer. Be aware that this date can be different from
        the "canonical" observation time (see the `.GenericMap.date` property).

        The reference date is determined using this order of preference:

        1. The ``DATE-AVG`` key in the meta.
        2. The ``DATE-OBS`` key in the meta.
        3. The ``DATE-BEG`` key in the meta.
        4. The ``DATE-END`` key in the meta.
        5. The `.GenericMap.date` property as a fallback (which, if not
           overridden, would be the current time if the above keywords are missing).

        See Also
        --------
        date : The observation time.
        date_start : The start time of the observation.
        date_end : The end time of the observation.
        date_average : The average time of the observation.

        Notes
        -----
        The FITS standard implies that, but does not expressly require, the DATE-AVG keyword
        to be the reference date.
        """
        return (
            self._get_date('date-avg') or
            self._date_obs or
            self._get_date('date-beg') or
            self._get_date('date-end') or
            self.date
        )

    def _set_reference_date(self, date):
        """
        Set the reference date using the same priority as `.GenericMap.reference_date`.

        If a source subclass overrides `.GenericMap.reference_date`, it should override
        this private method as well.
        """
        for keyword in ['date-avg', 'date-obs', 'date-beg', 'date-end']:
            if keyword in self.meta:
                self.meta[keyword] = parse_time(date).utc.isot
                return
        self._set_date(date)

    @property
    def date(self):
        """
        The observation time.

        This time is the "canonical" way to refer to an observation, which is commonly
        the start of the observation, but can be a different time. In comparison, the
        `.GenericMap.date_start` property is unambigiously the start of the observation.

        The observation time is determined using this order of preference:

        1. The ``DATE-OBS`` or ``DATE_OBS`` FITS keywords
        2. `.GenericMap.date_start`
        3. `.GenericMap.date_average`
        4. `.GenericMap.date_end`
        5. The current time

        See Also
        --------
        reference_date : The reference date for the the coordinate system
        date_start : The start time of the observation.
        date_end : The end time of the observation.
        date_average : The average time of the observation.
        """
        time = (
            self._date_obs or
            self.date_start or
            self.date_average or
            self.date_end
        )

        if time is None:
            if self._default_time is None:
                warn_metadata("Missing metadata for observation time, "
                              "setting observation time to current time. "
                              "Set the 'DATE-OBS' FITS keyword to prevent this warning.")
                self._default_time = parse_time('now')
            time = self._default_time

        return time

    def _set_date(self, date):
        """
        Set the observation time by setting DATE-OBS.

        If a source subclass overrides `.GenericMap.date`, it should override
        this private method as well.

        Notes
        -----
        This method will additionally always remove DATE_OBS (note the underscore),
        if present.
        """
        if 'date_obs' in self.meta:
            del self.meta['date_obs']
        self.meta['date-obs'] = parse_time(date).utc.isot

    @property
    def detector(self):
        """
        Detector name.

        This is taken from the 'DETECTOR' FITS keyword.
        """
        return self.meta.get('detector', "")

    @property
    def timeunit(self):
        """
        The `~astropy.units.Unit` of the exposure time of this observation.

        Taken from the "TIMEUNIT" FITS keyword, and defaults to seconds (as per)
        the FITS standard).
        """
        return u.Unit(self.meta.get('timeunit', 's'))

    @property
    def exposure_time(self):
        """
        Exposure time of the image.

        This is taken from the 'XPOSURE' keyword or the 'EXPTIME' FITS keyword,
        in that order.
        """
        exptime = self.meta.get('xposure') or self.meta.get('exptime')
        if exptime is not None:
            return exptime * self.timeunit

    @property
    def instrument(self):
        """Instrument name."""
        return self.meta.get('instrume', "").replace("_", " ")

    @property
    def measurement(self):
        """
        The measurement type of the observation.

        The measurement type can be described by a `str` or a
        `~astropy.units.Quantity`. If the latter, it is typically equal to
        `.GenericMap.wavelength`.

        See Also
        --------
        wavelength : The wavelength of the observation.
        """

        return self.wavelength

    @property
    def waveunit(self):
        """
        The `~astropy.units.Unit` of the wavelength of this observation.

        This is taken from the 'WAVEUNIT' FITS keyword. If the keyword is not
        present, defaults to `None`
        """
        if 'waveunit' in self.meta:
            return u.Unit(self.meta['waveunit'])
        else:
            wunit = extract_waveunit(self.meta)
            if wunit is not None:
                return u.Unit(wunit)

    @property
    def wavelength(self):
        """
        Wavelength of the observation.

        This is taken from the 'WAVELNTH' FITS keywords. If the keyword is not
        present, defaults to `None`. If 'WAVEUNIT' keyword isn't present,
        defaults to dimensionless units.
        """
        if 'wavelnth' in self.meta:
            return u.Quantity(self.meta['wavelnth'], self.waveunit)

    @property
    def observatory(self):
        """
        Observatory or Telescope name.

        This is taken from the 'OBSRVTRY' FITS keyword.
        """
        return self.meta.get('obsrvtry',
                             self.meta.get('telescop', "")).replace("_", " ")

    @property
    def processing_level(self):
        """
        Returns the FITS processing level if present.

        This is taken from the 'LVL_NUM' FITS keyword.
        """
        return self.meta.get('lvl_num', None)

    @property
    def bottom_left_coord(self):
        """
        The physical coordinate at the center of the bottom left ([0, 0]) pixel.
        """
        return self.wcs.pixel_to_world(0, 0)

    @property
    def top_right_coord(self):
        """
        The physical coordinate at the center of the the top right ([-1, -1]) pixel.
        """
        top_right = np.array([self.dimensions.x.value, self.dimensions.y.value]) - 1
        return self.wcs.pixel_to_world(*top_right)

    @property
    def center(self):
        """
        Return a coordinate object for the center pixel of the array.

        If the array has an even number of pixels in a given dimension,
        the coordinate returned lies on the edge between the two central pixels.
        """
        center = (np.array([self.dimensions.x.value, self.dimensions.y.value]) - 1) / 2.
        return self.wcs.pixel_to_world(*center)

    @u.quantity_input
    def shift_reference_coord(self, axis1: u.deg, axis2: u.deg):
        """
        Returns a map shifted by a specified amount to, for example, correct
        for a bad map location. These values are applied directly to the
        `~sunpy.map.GenericMap.reference_coordinate`. To check how much the
        reference coordinate has been modified, see
        ``sunpy.map.GenericMap.meta.modified_items['CRVAL1']`` and
        ``sunpy.map.GenericMap.meta.modified_items['CRVAL2']``.

        Parameters
        ----------
        axis1 : `~astropy.units.Quantity`
            The shift to apply to the Longitude (solar-x) coordinate.
        axis2 : `~astropy.units.Quantity`
            The shift to apply to the Latitude (solar-y) coordinate

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new shifted Map.
        """
        new_meta = self.meta.copy()

        # Update crvals
        new_meta['crval1'] = ((self._reference_longitude + axis1).to(self.spatial_units[0])).value
        new_meta['crval2'] = ((self._reference_latitude + axis2).to(self.spatial_units[1])).value

        # Create new map with the modification
        new_map = self._new_instance(self.data, new_meta, self.plot_settings)

        return new_map

    def _rsun_meters(self, dsun=None):
        """
        This property exists to avoid circular logic in constructing the
        observer coordinate, by allowing a custom 'dsun' to be specified,
        instead of one extracted from the `.observer_coordinate` property.
        """
        rsun = self.meta.get('rsun_ref', None)
        if rsun is not None:
            return rsun * u.m
        elif self._rsun_obs_no_default is not None:
            if dsun is None:
                dsun = self.dsun
            return sun._radius_from_angular_radius(self.rsun_obs, dsun)
        else:
            log.info("Missing metadata for solar radius: assuming "
                     "the standard radius of the photosphere.")
            return constants.radius

    @property
    def rsun_meters(self):
        """
        Assumed radius of observed emission from the Sun center.

        This is taken from the RSUN_REF FITS keyword, if present.
        If not, and angular radius metadata is present, it is calculated from
        `~sunpy.map.GenericMap.rsun_obs` and `~sunpy.map.GenericMap.dsun`.
        If neither pieces of metadata are present, defaults to the standard
        photospheric radius.
        """
        return self._rsun_meters()

    @property
    def _rsun_obs_no_default(self):
        """
        Get the angular radius value from FITS keywords without defaulting.
        Exists to avoid circular logic in `rsun_meters()` above.
        """
        return self.meta.get('rsun_obs',
                             self.meta.get('solar_r',
                                           self.meta.get('radius',
                                                         None)))

    @property
    def rsun_obs(self):
        """
        Angular radius of the observation from Sun center.

        This value is taken (in order of preference) from the 'RSUN_OBS',
        'SOLAR_R', or 'RADIUS' FITS keywords. If none of these keys are present,
        the angular radius is calculated from
        `~sunpy.map.GenericMap.rsun_meters` and `~sunpy.map.GenericMap.dsun`.
        """
        rsun_arcseconds = self._rsun_obs_no_default

        if rsun_arcseconds is not None:
            return rsun_arcseconds * u.arcsec
        else:
            return sun._angular_radius(self.rsun_meters, self.dsun)

    @property
    def coordinate_system(self):
        """
        Coordinate system used for x and y axes (ctype1/2).

        If not present, defaults to (HPLN-TAN, HPLT-TAN), and emits a warning.
        """
        ctype1 = self.meta.get('ctype1', None)
        if not ctype1:
            warn_metadata("Missing CTYPE1 from metadata, assuming CTYPE1 is HPLN-TAN")
            ctype1 = 'HPLN-TAN'

        ctype2 = self.meta.get('ctype2', None)
        if not ctype2:
            warn_metadata("Missing CTYPE2 from metadata, assuming CTYPE2 is HPLT-TAN")
            ctype2 = 'HPLT-TAN'

        # Astropy WCS does not understand the SOHO default of "solar-x" and
        # "solar-y" ctypes. This overrides the default assignment and
        # changes it to a ctype that is understood. See Thompson, 2006, A.&A.,
        # 449, 791.
        if ctype1.lower() in ("solar-x", "solar_x"):
            warn_deprecated("CTYPE1 value 'solar-x'/'solar_x' is deprecated, use 'HPLN-TAN' instead.")
            ctype1 = 'HPLN-TAN'

        if ctype2.lower() in ("solar-y", "solar_y"):
            warn_deprecated("CTYPE2 value 'solar-y'/'solar_y' is deprecated, use 'HPLN-TAN' instead.")
            ctype2 = 'HPLT-TAN'

        return SpatialPair(ctype1, ctype2)

    @property
    def _supported_observer_coordinates(self):
        """
        A list of supported coordinate systems.

        This is a list so it can easily maintain a strict order. The list of
        two element tuples, the first item in the tuple is the keys that need
        to be in the header to use this coordinate system and the second is the
        kwargs to SkyCoord.
        """
        return [(('hgln_obs', 'hglt_obs', 'dsun_obs'), {'lon': self.meta.get('hgln_obs'),
                                                        'lat': self.meta.get('hglt_obs'),
                                                        'radius': self.meta.get('dsun_obs'),
                                                        'unit': (u.deg, u.deg, u.m),
                                                        'frame': "heliographic_stonyhurst"}),
                (('crln_obs', 'crlt_obs', 'dsun_obs'), {'lon': self.meta.get('crln_obs'),
                                                        'lat': self.meta.get('crlt_obs'),
                                                        'radius': self.meta.get('dsun_obs'),
                                                        'unit': (u.deg, u.deg, u.m),
                                                        'frame': "heliographic_carrington"}), ]

    @property
    def _default_observer_coordinate(self):
        """
        The default observer coordinate to use when there is insufficient information
        in the metadata. This should be overridden by map sources as appropriate.
        """

    def _remove_existing_observer_location(self):
        """
        Remove all keys that this map might use for observer location.
        """
        all_keys = expand_list([e[0] for e in self._supported_observer_coordinates])
        for key in all_keys:
            self.meta.pop(key)

    @property
    @cached_property_based_on('_meta_hash')
    def observer_coordinate(self):
        """
        The Heliographic Stonyhurst Coordinate of the observer.

        Notes
        -----
        The ``obstime`` for this coordinate uses the `.reference_date` property, which
        may be different from the `.date` property.
        """
        warning_message = []
        for keys, kwargs in self._supported_observer_coordinates:
            missing_keys = set(keys) - self.meta.keys()
            if not missing_keys:
                sc = SkyCoord(obstime=self.reference_date, **kwargs)
                # If the observer location is supplied in Carrington coordinates,
                # the coordinate's `observer` attribute should be set to "self"
                if isinstance(sc.frame, HeliographicCarrington):
                    sc.frame._observer = "self"

                sc = sc.heliographic_stonyhurst
                # We set rsun after constructing the coordinate, as we need
                # the observer-Sun distance (sc.radius) to calculate this, which
                # may not be provided directly in metadata (if e.g. the
                # observer coordinate is specified in a cartesian
                # representation)
                return SkyCoord(sc.replicate(rsun=self._rsun_meters(sc.radius)))
            elif missing_keys != keys:
                frame = kwargs['frame'] if isinstance(kwargs['frame'], str) else kwargs['frame'].name
                warning_message.append(f"For frame '{frame}' the following metadata is missing: "
                                       f"{','.join(missing_keys)}")

        default = self._default_observer_coordinate
        if default is not None:
            # If a map source specifies a default observer, we log a message at the debug level
            warning_message = (["Missing metadata for observer: assuming custom default observer."]
                               + warning_message)
            log.debug("\n".join(warning_message))
            return default
        else:
            # If a map source does not specify a default observer, we assume Earth center and warn
            warning_message = (["Missing metadata for observer: assuming Earth-based observer."]
                               + warning_message + [""])
            warn_metadata("\n".join(warning_message), stacklevel=3)
            return get_earth(self.reference_date)

    @property
    def heliographic_latitude(self):
        """Observer heliographic latitude."""
        return self.observer_coordinate.lat

    @property
    def heliographic_longitude(self):
        """Observer heliographic longitude."""
        return self.observer_coordinate.lon

    @property
    def carrington_latitude(self):
        """Observer Carrington latitude."""
        hgc_frame = HeliographicCarrington(observer=self.observer_coordinate, obstime=self.reference_date,
                                           rsun=self.rsun_meters)
        return self.observer_coordinate.transform_to(hgc_frame).lat

    @property
    def carrington_longitude(self):
        """Observer Carrington longitude."""
        hgc_frame = HeliographicCarrington(observer=self.observer_coordinate, obstime=self.reference_date,
                                           rsun=self.rsun_meters)
        return self.observer_coordinate.transform_to(hgc_frame).lon

    @property
    def dsun(self):
        """Observer distance from the center of the Sun."""
        return self.observer_coordinate.radius.to('m')

    @property
    def _reference_longitude(self):
        """
        FITS-WCS compatible longitude. Used in self.wcs and
        self.reference_coordinate.
        """
        return self.meta.get('crval1', 0.) * self.spatial_units[0]

    @property
    def _reference_latitude(self):
        return self.meta.get('crval2', 0.) * self.spatial_units[1]

    @property
    def reference_coordinate(self):
        """Reference point WCS axes in data units (i.e. crval1, crval2). This value
        includes a shift if one is set."""
        return SkyCoord(self._reference_longitude,
                        self._reference_latitude,
                        frame=self.coordinate_frame)

    @property
    def reference_pixel(self):
        """
        Pixel of reference coordinate.

        The pixel returned uses zero-based indexing, so will be 1 pixel less
        than the FITS CRPIX values.
        """
        naxis1 = self.meta.get('naxis1', self.data.shape[1])
        naxis2 = self.meta.get('naxis2', self.data.shape[0])
        return PixelPair((self.meta.get('crpix1', (naxis1 + 1) / 2.) - 1) * u.pixel,
                         (self.meta.get('crpix2', (naxis2 + 1) / 2.) - 1) * u.pixel)

    @property
    def scale(self):
        """
        Image scale along the x and y axes in units/pixel
        (i.e. cdelt1, cdelt2).

        If the CDij matrix is defined but no CDELTi values are explicitly defined,
        effective CDELTi values are constructed from the CDij matrix. The effective
        CDELTi values are chosen so that each row of the PCij matrix has unity norm.
        This choice is optimal if the PCij matrix is a pure rotation matrix, but may not
        be as optimal if the PCij matrix includes any skew.
        """
        if 'cd1_1' in self.meta and 'cdelt1' not in self.meta and 'cdelt2' not in self.meta:
            cdelt1 = np.sqrt(self.meta['cd1_1']**2 + self.meta['cd1_2']**2)
            cdelt2 = np.sqrt(self.meta['cd2_1']**2 + self.meta['cd2_2']**2)
        else:
            cdelt1 = self.meta.get('cdelt1', 1.)
            cdelt2 = self.meta.get('cdelt2', 1.)

        return SpatialPair(cdelt1 * self.spatial_units[0] / u.pixel,
                           cdelt2 * self.spatial_units[1] / u.pixel)

    @property
    def spatial_units(self):
        """
        Image coordinate units along the x and y axes (i.e. cunit1, cunit2).
        """
        units = self.meta.get('cunit1', None), self.meta.get('cunit2', None)
        units = [None if unit is None else u.Unit(unit.lower()) for unit in units]
        return SpatialPair(units[0], units[1])

    @property
    def rotation_matrix(self):
        r"""
        Matrix describing the transformation needed to align the reference
        pixel with the coordinate axes.

        The order or precedence of FITS keywords which this is taken from is:
        - PC\*_\*
        - CD\*_\*
        - CROTA\*

        Notes
        -----
        In many cases this is a simple rotation matrix, hence the property name.
        It general it does not have to be a pure rotation matrix, and can encode
        other transformations e.g., skews for non-orthogonal coordinate systems.
        """
        if any(key in self.meta for key in ['PC1_1', 'PC1_2', 'PC2_1', 'PC2_2']):
            return np.array(
                [
                    [self.meta.get('PC1_1', 1), self.meta.get('PC1_2', 0)],
                    [self.meta.get('PC2_1', 0), self.meta.get('PC2_2', 1)]
                ]
            )
        elif any(key in self.meta for key in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']):
            cd = np.array(
                [
                    [self.meta.get('CD1_1', 0), self.meta.get('CD1_2', 0)],
                    [self.meta.get('CD2_1', 0), self.meta.get('CD2_2', 0)]
                ]
            )

            cdelt = u.Quantity(self.scale).value

            # Divide each row by each CDELT
            return cd / np.expand_dims(cdelt, axis=1)
        else:
            return self._rotation_matrix_from_crota()

    @staticmethod
    def _pc_matrix(lam, angle):
        """
        Returns PC matrix from the scale ration (lam) and rotation
        angle in radians (angle).
        """
        return np.array([[np.cos(angle), -1 * lam * np.sin(angle)],
                         [1/lam * np.sin(angle), np.cos(angle)]])

    def _rotation_matrix_from_crota(self, crota_key='CROTA2'):
        """
        This method converts the deprecated CROTA FITS kwargs to the new
        PC rotation matrix.

        This method can be overridden if an instruments header does not use this
        conversion.

        Parameters
        ----------
        crota_key : str, optional
            The key to use for CROTA2. Defaults to 'CROTA2'.

        Notes
        -----
        If the specified key isn't present in the metadata, a default rotation
        of 0deg is returned.
        """
        lam = self.scale[1] / self.scale[0]
        p = np.deg2rad(self.meta.get(crota_key, 0))
        return self._pc_matrix(lam, p)

    @property
    def _pv_values(self):
        """
        Return any PV values in the metadata.
        """
        pattern = re.compile(r'pv[1-9]\d?_(?:0|[1-9]\d?)$', re.IGNORECASE)
        pv_keys = [k for k in self.meta.keys() if pattern.match(k)]

        pv_values = []
        for k in pv_keys:
            i, m = int(k[2]), int(k[4:])
            pv_values.append((i, m, self.meta[k]))
        return pv_values

    @property
    def fits_header(self):
        """
        A `~astropy.io.fits.Header` representation of the ``meta`` attribute.
        """
        return header_to_fits(self.meta)

# #### Miscellaneous #### #
    def _get_cmap_name(self):
        """Build the default color map name."""
        cmap_string = (self.observatory + self.detector +
                       str(int(self.wavelength.to('angstrom').value)))
        return cmap_string.lower()

    def _validate_meta(self):
        """
        Validates some meta-information associated with a Map.

        This method includes very basic validation checks which apply to
        all of the kinds of files that sunpy can read. Datasource-specific
        validation should be handled in the relevant file in the
        sunpy.map.sources package.
        """
        msg = ('Image coordinate units for axis {} not present in metadata.')
        err_message = []
        for i in [0, 1]:
            if self.spatial_units[i] is None:
                err_message.append(msg.format(i+1, i+1))

        if err_message:
            err_message.append(
                f'See {_META_FIX_URL} for instructions on how to add missing metadata.')
            raise MapMetaValidationError('\n'.join(err_message))

        for meta_property in ('waveunit', ):
            if (self.meta.get(meta_property) and
                u.Unit(self.meta.get(meta_property),
                       parse_strict='silent').physical_type == 'unknown'):
                warn_metadata(f"Unknown value for {meta_property.upper()}.")

        if (self.coordinate_system[0].startswith(('SOLX', 'SOLY')) or
                self.coordinate_system[1].startswith(('SOLX', 'SOLY'))):
            warn_user("sunpy Map does not support three dimensional data "
                      "and therefore cannot represent heliocentric coordinates. Proceed at your own risk.")

        if not all(su.is_equivalent(u.arcsec) for su in self.spatial_units):
            units = [su.to_string() for su in self.spatial_units]
            raise MapMetaValidationError(
                'Map only supports spherical coordinate systems with angular units '
                f'(ie. equivalent to arcsec), but this map has units {units}')

# #### Data conversion routines #### #
    def world_to_pixel(self, coordinate):
        """
        Convert a world (data) coordinate to a pixel coordinate.

        Parameters
        ----------
        coordinate : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseCoordinateFrame`
            The coordinate object to convert to pixel coordinates.

        Returns
        -------
        x : `~astropy.units.Quantity`
            Pixel coordinate on the CTYPE1 axis.
        y : `~astropy.units.Quantity`
            Pixel coordinate on the CTYPE2 axis.
        """
        x, y = self.wcs.world_to_pixel(coordinate)
        return PixelPair(x * u.pixel, y * u.pixel)

    @u.quantity_input
    def pixel_to_world(self, x: u.pixel, y: u.pixel):
        """
        Convert a pixel coordinate to a data (world) coordinate.

        Parameters
        ----------
        x : `~astropy.units.Quantity`
            Pixel coordinate of the CTYPE1 axis. (Normally solar-x).
        y : `~astropy.units.Quantity`
            Pixel coordinate of the CTYPE2 axis. (Normally solar-y).

        Returns
        -------
        coord : `astropy.coordinates.SkyCoord`
            A coordinate object representing the output coordinate.
        """
        return self.wcs.pixel_to_world(x, y)

# #### I/O routines #### #

    def save(self, filepath, filetype='auto', **kwargs):
        """
        Save a map to a file.

        Parameters
        ----------
        filepath : `str`
            Location to save the file to.
            If ``filepath`` ends with ".asdf" and ``filetype="auto"``, an ASDF file will be created.
        filetype : `str`, optional
            The file format to save the map in. Defaults to ``"auto"`` which infers
            the format from the file extension. Supported formats include FITS, JP2, and ASDF.
        hdu_type : `~astropy.io.fits.hdu.base.ExtensionHDU` instance or class, optional
            For FITS files, this specifies the type of HDU to write. By default, the map is saved
            in the primary HDU. If an HDU type or instance is provided, the map data and header will
            be written to that HDU. For example, `astropy.io.fits.CompImageHDU` can be used to compress the map.
        kwargs :
            Any additional keyword arguments are passed to `~sunpy.io._file_tools.write_file`
            or `asdf.AsdfFile.write_to`.

        Notes
        -----
        Saving with the jp2 extension will write a modified version
        of the given data casted to uint8 values in order to support
        the JPEG2000 format.

        Saving with the ``.asdf`` extension will save the map as an ASDF file, storing the map's
        attributes under the key ``'sunpymap'`` in the ASDF tree.

        Examples
        --------
        >>> from astropy.io.fits import CompImageHDU
        >>> from sunpy.map import Map
        >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
        >>> aia_map = Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
        >>> aia_map.save("aia171.fits", hdu_type=CompImageHDU)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
        >>> aia_map.save("aia171.asdf")  # doctest: +REMOTE_DATA
        """
        if filetype.lower() == "asdf" or (filetype.lower() == "auto" and str(filepath).lower().endswith(".asdf")):
            import asdf
            asdf.AsdfFile({'sunpymap': self}).write_to(str(filepath), **kwargs)
        else:
            write_file(filepath, self.data, self.meta, filetype=filetype, **kwargs)




# #### Image processing routines #### #

    @u.quantity_input
    def resample(self, dimensions: u.pixel, method='linear'):
        """
        Resample to new dimension sizes.

        Uses the same parameters and creates the same coordinate lookup points
        as IDL''s congrid routine, which apparently originally came from a
        VAX/VMS routine of the same name.

        This method **does not** preserve dask arrays.

        Parameters
        ----------
        dimensions : `~astropy.units.Quantity`
            Output pixel dimensions. The first argument corresponds to the 'x'
            axis and the second argument corresponds to the 'y' axis.
        method : str
            Method to use for resampling interpolation.
                * ``'nearest'`` and ``'linear'`` - Use n x 1-D interpolations using
                  `scipy.interpolate.interp1d`.
                * ``'spline'`` - Uses piecewise polynomials (splines) for mapping the input
                  array to new coordinates by interpolation using
                  `scipy.ndimage.map_coordinates`.

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            Resampled map

        References
        ----------
        `Rebinning <https://scipy-cookbook.readthedocs.io/items/Rebinning.html>`__
        """

        # Note: because the underlying ndarray is transposed in sense when
        #   compared to the Map, the ndarray is transposed, resampled, then
        #   transposed back
        # Note: "center" defaults to True in this function because data
        #   coordinates in a Map are at pixel centers

        # Make a copy of the original data and perform resample
        new_data = sunpy_image_resample(self.data.copy().T, dimensions,
                                        method, center=True)
        new_data = new_data.T

        scale_factor_x = float(self.dimensions[0] / dimensions[0])
        scale_factor_y = float(self.dimensions[1] / dimensions[1])

        # Update image scale and number of pixels
        new_meta = self.meta.copy()

        # Update metadata
        for key in {'cdelt1', 'cd1_1', 'cd2_1'} & self.meta.keys():
            new_meta[key] *= scale_factor_x
        for key in {'cdelt2', 'cd1_2', 'cd2_2'} & self.meta.keys():
            new_meta[key] *= scale_factor_y
        if 'pc1_1' in self.meta:
            new_meta['pc1_2'] *= scale_factor_y / scale_factor_x
            new_meta['pc2_1'] *= scale_factor_x / scale_factor_y

        new_meta['crpix1'] = (self.reference_pixel.x.to_value(u.pix) + 0.5) / scale_factor_x + 0.5
        new_meta['crpix2'] = (self.reference_pixel.y.to_value(u.pix) + 0.5) / scale_factor_y + 0.5
        new_meta['naxis1'] = new_data.shape[1]
        new_meta['naxis2'] = new_data.shape[0]

        # Create new map instance
        new_map = self._new_instance(new_data, new_meta, self.plot_settings)
        return new_map

    @add_common_docstring(rotation_function_names=_rotation_function_names)
    @u.quantity_input
    def rotate(self, angle: u.deg = None, rmatrix=None, order=3, scale=1.0,
               recenter=False, missing=np.nan, *, method='scipy', clip=True):
        """
        Returns a new rotated and rescaled map.

        Specify either a rotation angle or a rotation matrix, but not both. If
        neither an angle or a rotation matrix are specified, the map will be
        rotated by the rotation information in the metadata, which should derotate
        the map so that the pixel axes are aligned with world-coordinate axes.

        This method **does not** preserve dask arrays.

        Parameters
        ----------
        angle : `~astropy.units.Quantity`
            The angle (degrees) to rotate counterclockwise.
        rmatrix : array-like
            2x2 linear transformation rotation matrix.
        order : int
            Interpolation order to be used. The precise meaning depends on the
            rotation method specified by ``method``.
            Default: 3
        scale : float
            A scale factor for the image, default is no scaling
        recenter : bool
            If `True`, position the reference coordinate at the center of the new map
            Default: `False`
        missing : float
            The value to use for pixels in the output map that are beyond the extent
            of the input map.
            Default: `numpy.nan`
        method : {{{rotation_function_names}}}, optional
            Rotation function to use. Defaults to ``'scipy'``.
        clip : `bool`, optional
            If `True`, clips the pixel values of the output image to the range of the
            input image (including the value of ``missing``, if used).
            Defaults to `True`.

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new Map instance containing the rotated and rescaled data of the
            original map.

        See Also
        --------
        sunpy.image.transform.affine_transform :
            The routine this method calls for the rotation.

        Notes
        -----
        The rotation information in the metadata of the new map is modified
        appropriately from that of the original map to account for the applied rotation.
        It will solely be in the form of a PCi_j matrix, even if the original map used
        the CROTA2 keyword or a CDi_j matrix.

        If the map does not already contain pixels with `numpy.nan` values, setting
        ``missing`` to an appropriate number for the data (e.g., zero) will reduce the
        computation time.

        For each NaN pixel in the input image, one or more pixels in the output image
        will be set to NaN, with the size of the pixel region affected depending on the
        interpolation order. All currently implemented rotation methods require a
        convolution step to handle image NaNs. This convolution normally uses
        :func:`scipy.signal.convolve2d`, but if `OpenCV <https://opencv.org>`__ is
        installed, the faster |cv2_filter2D|_ is used instead.

        See :func:`sunpy.image.transform.affine_transform` for details on each of the
        rotation functions.
        """
        if angle is not None and rmatrix is not None:
            raise ValueError("You cannot specify both an angle and a rotation matrix.")
        elif angle is None and rmatrix is None:
            # Be aware that self.rotation_matrix may not actually be a pure rotation matrix
            rmatrix = self.rotation_matrix

        if order not in range(6):
            raise ValueError("Order must be between 0 and 5.")

        method = _get_transform_method(method)

        # The FITS-WCS transform is by definition defined around the
        # reference coordinate in the header.
        lon, lat = self._get_lon_lat(self.reference_coordinate.frame)
        rotation_center = u.Quantity([lon, lat])

        # Copy meta data
        new_meta = self.meta.copy()
        if angle is not None:
            # Calculate the parameters for the affine_transform
            c = np.cos(np.deg2rad(angle))
            s = np.sin(np.deg2rad(angle))
            rmatrix = np.array([[c, -s],
                                [s, c]])

        # The data will be rotated by the inverse of the rotation matrix. Because rmatrix may not
        # actually be a pure rotation matrix, we calculate the inverse in a general manner.
        inv_rmatrix = np.linalg.inv(rmatrix)

        # Calculate the shape in pixels to contain all of the image data
        corners = itertools.product([-0.5, self.data.shape[1]-0.5],
                                    [-0.5, self.data.shape[0]-0.5])
        rot_corners = np.vstack([rmatrix @ c for c in corners]) * scale
        extent = np.max(rot_corners, axis=0) - np.min(rot_corners, axis=0)

        # Calculate the needed padding or unpadding
        diff = np.asarray(np.ceil((extent - np.flipud(self.data.shape)) / 2), dtype=int)
        pad_x = np.max((diff[0], 0))
        pad_y = np.max((diff[1], 0))
        unpad_x = -np.min((diff[0], 0))
        unpad_y = -np.min((diff[1], 0))

        # Raise an informative error message if trying to pad an integer array with NaNs
        if (pad_x > 0 or pad_y > 0) and issubclass(self.data.dtype.type, numbers.Integral) and (missing % 1 != 0):
            raise ValueError("The underlying data is integers, but the fill value for missing "
                             "pixels cannot be cast to an integer, which is the case for the "
                             "default fill value of NaN. Set the `missing` keyword to an "
                             "appropriate integer value for the data set.")

        new_data = np.pad(self.data,
                          ((pad_y, pad_y), (pad_x, pad_x)),
                          mode='constant',
                          constant_values=(missing, missing))

        # All of the following pixel calculations use a pixel origin of 0
        pixel_array_center = (np.flipud(new_data.shape) - 1) / 2.0
        pixel_rotation_center = u.Quantity(self.reference_pixel).value + [pad_x, pad_y]

        if recenter:
            pixel_center = pixel_rotation_center
        else:
            pixel_center = pixel_array_center

        # Apply the rotation to the image data
        new_data = affine_transform(new_data,
                                    np.asarray(inv_rmatrix),
                                    order=order, scale=scale,
                                    image_center=pixel_center,
                                    recenter=recenter, missing=missing,
                                    method=method, clip=clip)

        if recenter:
            new_reference_pixel = pixel_array_center
        else:
            # Calculate new pixel coordinates for the rotation center
            new_reference_pixel = pixel_center + np.dot(rmatrix * scale,
                                                        pixel_rotation_center - pixel_center)
            new_reference_pixel = np.array(new_reference_pixel).ravel()

        # Define the new reference_pixel
        new_meta['crval1'] = rotation_center[0].value
        new_meta['crval2'] = rotation_center[1].value
        new_meta['crpix1'] = new_reference_pixel[0] + 1  # FITS pixel origin is 1
        new_meta['crpix2'] = new_reference_pixel[1] + 1  # FITS pixel origin is 1
        new_meta['NAXIS1'] = new_data.shape[1]
        new_meta['NAXIS2'] = new_data.shape[0]

        # Unpad the array if necessary
        if unpad_x > 0:
            new_data = new_data[:, unpad_x:-unpad_x]
            new_meta['crpix1'] -= unpad_x
        if unpad_y > 0:
            new_data = new_data[unpad_y:-unpad_y, :]
            new_meta['crpix2'] -= unpad_y

        # Calculate the new rotation matrix to store in the header by
        # "subtracting" the rotation matrix used in the rotate from the old one
        # That being calculate the dot product of the old header data with the
        # inverse of the rotation matrix.
        pc_C = np.dot(self.rotation_matrix, inv_rmatrix)
        new_meta['PC1_1'] = pc_C[0, 0]
        new_meta['PC1_2'] = pc_C[0, 1]
        new_meta['PC2_1'] = pc_C[1, 0]
        new_meta['PC2_2'] = pc_C[1, 1]

        # Update pixel size if image has been scaled.
        if scale != 1.0:
            new_meta['cdelt1'] = (self.scale[0] / scale).value
            new_meta['cdelt2'] = (self.scale[1] / scale).value

        # Remove old CROTA kwargs because we have saved a new PCi_j matrix.
        new_meta.pop('CROTA1', None)
        new_meta.pop('CROTA2', None)
        # Remove CDi_j header
        new_meta.pop('CD1_1', None)
        new_meta.pop('CD1_2', None)
        new_meta.pop('CD2_1', None)
        new_meta.pop('CD2_2', None)

        # Create new map with the modification
        new_map = self._new_instance(new_data, new_meta, self.plot_settings)

        return new_map

    @u.quantity_input
    def submap(self, bottom_left, *, top_right=None, width: (u.deg, u.pix) = None, height: (u.deg, u.pix) = None):
        """
        Returns a submap defined by a rectangle.

        Any pixels which have at least part of their area inside the rectangle
        are returned. If the rectangle is defined in world coordinates, the
        smallest array which contains all four corners of the rectangle as
        defined in world coordinates is returned.

        This method **does** preserve dask arrays.

        Parameters
        ----------
        bottom_left : `astropy.units.Quantity` or `~astropy.coordinates.SkyCoord`
            The bottom-left coordinate of the rectangle. If a `~astropy.coordinates.SkyCoord` it can
            have shape ``(2,)`` and simultaneously define ``top_right``. If specifying
            pixel coordinates it must be given as an `~astropy.units.Quantity`
            object with units of pixels.
        top_right : `astropy.units.Quantity` or `~astropy.coordinates.SkyCoord`, optional
            The top-right coordinate of the rectangle. If ``top_right`` is
            specified ``width`` and ``height`` must be omitted.
        width : `astropy.units.Quantity`, optional
            The width of the rectangle. Required if ``top_right`` is omitted.
        height : `astropy.units.Quantity`
            The height of the rectangle. Required if ``top_right`` is omitted.

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new map instance is returned representing to specified
            sub-region.

        Notes
        -----
        When specifying pixel coordinates, they are specified in Cartesian
        order not in numpy order. So, for example, the ``bottom_left=``
        argument should be ``[left, bottom]``.

        Examples
        --------
        >>> import astropy.units as u
        >>> from astropy.coordinates import SkyCoord
        >>> import sunpy.map
        >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
        >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
        >>> bl = SkyCoord(-300*u.arcsec, -300*u.arcsec, frame=aia.coordinate_frame)  # doctest: +REMOTE_DATA
        >>> tr = SkyCoord(500*u.arcsec, 500*u.arcsec, frame=aia.coordinate_frame)  # doctest: +REMOTE_DATA
        >>> aia.submap(bl, top_right=tr)   # doctest: +REMOTE_DATA
        <sunpy.map.sources.sdo.AIAMap object at ...>
        SunPy Map
        ---------
        Observatory:         SDO
        Instrument:          AIA 3
        Detector:            AIA
        Measurement:         171.0 Angstrom
        Wavelength:          171.0 Angstrom
        Observation Date:    2011-06-07 06:33:02
        Exposure Time:       0.234256 s
        Dimension:           [335. 335.] pix
        Coordinate System:   helioprojective
        Scale:               [2.402792 2.402792] arcsec / pix
        Reference Pixel:     [126.5 125.5] pix
        Reference Coord:     [3.22309951 1.38578135] arcsec
        ...

        >>> aia.submap([0,0]*u.pixel, top_right=[5,5]*u.pixel)   # doctest: +REMOTE_DATA
        <sunpy.map.sources.sdo.AIAMap object at ...>
        SunPy Map
        ---------
        Observatory:         SDO
        Instrument:          AIA 3
        Detector:            AIA
        Measurement:         171.0 Angstrom
        Wavelength:          171.0 Angstrom
        Observation Date:    2011-06-07 06:33:02
        Exposure Time:       0.234256 s
        Dimension:           [6. 6.] pix
        Coordinate System:   helioprojective
        Scale:               [2.402792 2.402792] arcsec / pix
        Reference Pixel:     [511.5 511.5] pix
        Reference Coord:     [3.22309951 1.38578135] arcsec
        ...

        >>> width = 10 * u.arcsec
        >>> height = 10 * u.arcsec
        >>> aia.submap(bl, width=width, height=height)   # doctest: +REMOTE_DATA
        <sunpy.map.sources.sdo.AIAMap object at ...>
        SunPy Map
        ---------
        Observatory:         SDO
        Instrument:          AIA 3
        Detector:            AIA
        Measurement:         171.0 Angstrom
        Wavelength:          171.0 Angstrom
        Observation Date:    2011-06-07 06:33:02
        Exposure Time:       0.234256 s
        Dimension:           [5. 5.] pix
        Coordinate System:   helioprojective
        Scale:               [2.402792 2.402792] arcsec / pix
        Reference Pixel:     [125.5 125.5] pix
        Reference Coord:     [3.22309951 1.38578135] arcsec
        ...

        >>> bottom_left_vector = SkyCoord([0, 10]  * u.deg, [0, 10] * u.deg, frame='heliographic_stonyhurst')
        >>> aia.submap(bottom_left_vector)   # doctest: +REMOTE_DATA
        <sunpy.map.sources.sdo.AIAMap object at ...>
        SunPy Map
        ---------
        Observatory:         SDO
        Instrument:          AIA 3
        Detector:            AIA
        Measurement:         171.0 Angstrom
        Wavelength:          171.0 Angstrom
        Observation Date:    2011-06-07 06:33:02
        Exposure Time:       0.234256 s
        Dimension:           [70. 69.] pix
        Coordinate System:   helioprojective
        Scale:               [2.402792 2.402792] arcsec / pix
        Reference Pixel:     [1.5 0.5] pix
        Reference Coord:     [3.22309951 1.38578135] arcsec
        ...
        """
        # Check that we have been given a valid combination of inputs
        # [False, False, False] is valid if bottom_left contains the two corner coords
        if ([arg is not None for arg in (top_right, width, height)]
                not in [[True, False, False], [False, False, False], [False, True, True]]):
            raise ValueError("Either top_right alone or both width and height must be specified.")
        # parse input arguments
        pixel_corners = u.Quantity(self._parse_submap_input(
            bottom_left, top_right, width, height)).T

        msg = (
            "The provided input coordinates to ``submap`` when transformed to the target "
            "coordinate frame contain NaN values and cannot be used to crop the map. "
            "The most common reason for NaN values is transforming off-disk 2D "
            "coordinates without specifying an assumption (e.g., via the"
            "`sunpy.coordinates.SphericalScreen()` context manager) that allows "
            "such coordinates to be interpreted as 3D coordinates."
        )
        if np.any(np.isnan(pixel_corners)):
            raise ValueError(msg)

        # The pixel corners result is in Cartesian order, so the first index is
        # columns and the second is rows.
        bottom = np.min(pixel_corners[1]).to_value(u.pix)
        top = np.max(pixel_corners[1]).to_value(u.pix)
        left = np.min(pixel_corners[0]).to_value(u.pix)
        right = np.max(pixel_corners[0]).to_value(u.pix)

        # Round the lower left pixel to the nearest integer
        # We want 0.5 to be rounded up to 1, so use floor(x + 0.5)
        bottom = np.floor(bottom + 0.5)
        left = np.floor(left + 0.5)

        # Round the top right pixel to the nearest integer, then add 1 for array indexing
        # We want e.g. 2.5 to be rounded down to 2, so use ceil(x - 0.5)
        top = np.ceil(top - 0.5) + 1
        right = np.ceil(right - 0.5) + 1

        # Clip pixel values to max of array, prevents negative
        # indexing
        bottom = int(np.clip(bottom, 0, self.data.shape[0]))
        top = int(np.clip(top, 0, self.data.shape[0]))
        left = int(np.clip(left, 0, self.data.shape[1]))
        right = int(np.clip(right, 0, self.data.shape[1]))

        arr_slice = np.s_[bottom:top, left:right]
        # Get ndarray representation of submap
        new_data = self.data[arr_slice].copy()

        # Make a copy of the header with updated centering information
        new_meta = self.meta.copy()
        # Add one to go from zero-based to one-based indexing
        new_meta['crpix1'] = self.reference_pixel.x.to_value(u.pix) + 1 - left
        new_meta['crpix2'] = self.reference_pixel.y.to_value(u.pix) + 1 - bottom
        new_meta['naxis1'] = new_data.shape[1]
        new_meta['naxis2'] = new_data.shape[0]

        # Create new map instance
        if self.mask is not None:
            new_mask = self.mask[arr_slice].copy()
            # Create new map with the modification
            new_map = self._new_instance(new_data, new_meta, self.plot_settings, mask=new_mask)
            return new_map
        # Create new map with the modification
        new_map = self._new_instance(new_data, new_meta, self.plot_settings)
        return new_map

    @seconddispatch
    def _parse_submap_input(self, bottom_left, top_right, width, height):
        """
        Should take any valid input to submap() and return bottom_left and
        top_right in pixel coordinates.
        """

    @_parse_submap_input.register(u.Quantity)
    def _parse_submap_quantity_input(self, bottom_left, top_right, width, height):
        if top_right is None and width is None:
            raise ValueError('Either top_right alone or both width and height must be specified '
                             'when bottom_left is a Quantity')
        if bottom_left.shape != (2, ):
            raise ValueError('bottom_left must have shape (2, ) when specified as a Quantity')
        if top_right is not None:
            if top_right.shape != (2, ):
                raise ValueError('top_right must have shape (2, ) when specified as a Quantity')
            if not top_right.unit.is_equivalent(u.pix):
                raise TypeError("When bottom_left is a Quantity, top_right "
                                "must be a Quantity in units of pixels.")
            # Have bottom_left and top_right in pixels already, so no need to do
            # anything else
        else:
            if not (width.unit.is_equivalent(u.pix) and
                    height.unit.is_equivalent(u.pix)):
                raise TypeError("When bottom_left is a Quantity, width and height "
                                "must be a Quantity in units of pixels.")
            # Add width and height to get top_right
            top_right = u.Quantity([bottom_left[0] + width, bottom_left[1] + height])

        top_left = u.Quantity([top_right[0], bottom_left[1]])
        bottom_right = u.Quantity([bottom_left[0], top_right[1]])
        return bottom_left, top_left, top_right, bottom_right

    @_parse_submap_input.register(SkyCoord)
    @_parse_submap_input.register(BaseCoordinateFrame)
    def _parse_submap_coord_input(self, bottom_left, top_right, width, height):
        # Use helper function to get top_right as a SkyCoord
        bottom_left, top_right = get_rectangle_coordinates(bottom_left,
                                                           top_right=top_right,
                                                           width=width,
                                                           height=height)

        if isinstance(bottom_left, SkyCoord):
            frame = bottom_left.frame

        frame = bottom_left
        left_lon, bottom_lat = self._get_lon_lat(bottom_left)
        right_lon, top_lat = self._get_lon_lat(top_right)
        corners = SkyCoord([left_lon, left_lon, right_lon, right_lon],
                           [bottom_lat, top_lat, top_lat, bottom_lat],
                           frame=frame)

        return tuple(u.Quantity(self.wcs.world_to_pixel(corners), u.pix).T)

    @u.quantity_input
    def superpixel(self, dimensions: u.pixel, offset: u.pixel = (0, 0)*u.pixel, func=np.sum,
                   conservative_mask: bool = False):
        """Returns a new map consisting of superpixels formed by applying
        'func' to the original map data.

        This method **does** preserve dask arrays.

        Parameters
        ----------
        dimensions : tuple
            One superpixel in the new map is equal to (dimension[0],
            dimension[1]) pixels of the original map.
            The first argument corresponds to the 'x' axis and the second
            argument corresponds to the 'y' axis. If non-integer values are provided,
            they are rounded using `int`.
        offset : tuple
            Offset from (0,0) in original map pixels used to calculate where
            the data used to make the resulting superpixel map starts.
            If non-integer value are provided, they are rounded using `int`.
        func
            Function applied to the original data.
            The function 'func' must take a numpy array as its first argument,
            and support the axis keyword with the meaning of a numpy axis
            keyword (see the description of `~numpy.sum` for an example.)
            The default value of 'func' is `~numpy.sum`; using this causes
            superpixel to sum over (dimension[0], dimension[1]) pixels of the
            original map.
        conservative_mask : bool, optional
            If `True`, a superpixel is masked if any of its constituent pixels are masked.
            If `False`, a superpixel is masked only if all of its constituent pixels are masked.
            Default is `False`.

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new Map which has superpixels of the required size.
        """

        # Note: because the underlying ndarray is transposed in sense when
        #   compared to the Map, the ndarray is transposed, resampled, then
        #   transposed back.
        # Note: "center" defaults to True in this function because data
        #   coordinates in a Map are at pixel centers.

        if (offset.value[0] < 0) or (offset.value[1] < 0):
            raise ValueError("Offset is strictly non-negative.")

        # These are rounded by int() in reshape_image_to_4d_superpixel,
        # so round here too for use in constructing metadata later.
        dimensions = [int(dim) for dim in dimensions.to_value(u.pix)]
        offset = [int(off) for off in offset.to_value(u.pix)]

        # Make a copy of the original data, perform reshaping, and apply the
        # function.
        if self.mask is not None:
            data = np.ma.array(self.data.copy(), mask=self.mask)
        else:
            data = self.data.copy()

        reshaped_data = reshape_image_to_4d_superpixel(data, [dimensions[1], dimensions[0]], [offset[1], offset[0]])
        new_array = func(func(reshaped_data, axis=3), axis=1)

        if self.mask is not None:
            if conservative_mask ^ (func in [np.sum, np.prod]):
                log.info(
                    f"Using conservative_mask={conservative_mask} for function {func.__name__}, "
                    "which may not be ideal. Recommended: conservative_mask=True for sum/prod, "
                    "False for mean/median/std/min/max."
                    )

            if conservative_mask:
                reshaped_mask = reshape_image_to_4d_superpixel(self.mask, [dimensions[1], dimensions[0]], [offset[1], offset[0]])
                new_mask = np.any(reshaped_mask, axis=(3, 1))
            else:
                new_mask = np.ma.getmaskarray(new_array)
        else:
            new_mask = None

        # Update image scale and number of pixels

        # create copy of new meta data
        new_meta = self.meta.copy()

        # Update metadata
        for key in {'cdelt1', 'cd1_1', 'cd2_1'} & self.meta.keys():
            new_meta[key] *= dimensions[0]
        for key in {'cdelt2', 'cd1_2', 'cd2_2'} & self.meta.keys():
            new_meta[key] *= dimensions[1]
        if 'pc1_1' in self.meta:
            new_meta['pc1_2'] *= dimensions[1] / dimensions[0]
            new_meta['pc2_1'] *= dimensions[0] / dimensions[1]

        new_meta['crpix1'] = ((self.reference_pixel.x.to_value(u.pix) +
                               0.5 - offset[0]) / dimensions[0]) + 0.5
        new_meta['crpix2'] = ((self.reference_pixel.y.to_value(u.pix) +
                               0.5 - offset[1]) / dimensions[1]) + 0.5
        new_meta['naxis1'] = new_array.shape[1]
        new_meta['naxis2'] = new_array.shape[0]

        # Create new map instance
        if self.mask is not None:
            new_data = np.ma.getdata(new_array)
        else:
            new_data = new_array

        # Create new map with the modified data
        new_map = self._new_instance(new_data, new_meta, self.plot_settings, mask=new_mask)
        return new_map

# #### Visualization #### #

    @property
    def cmap(self):
        """
        Return the `matplotlib.colors.Colormap` instance this map uses.
        """
        cmap = self.plot_settings['cmap']
        if isinstance(cmap, str):
            cmap = plt.get_cmap(cmap)
            # Set the colormap to be this specific instance so we are not
            # returning a copy
            self.plot_settings['cmap'] = cmap
        return cmap

    @u.quantity_input
    def draw_grid(self, axes=None, grid_spacing: u.deg = 15*u.deg, annotate=True,
                  system='stonyhurst', **kwargs):
        """
        Draws a coordinate overlay on the plot in the Heliographic Stonyhurst
        coordinate system.

        To overlay other coordinate systems see the `WCSAxes Documentation
        <https://docs.astropy.org/en/stable/visualization/wcsaxes/overlaying_coordinate_systems.html>`__

        Parameters
        ----------
        axes : `~matplotlib.axes` or `None`
            Axes to plot limb on, or `None` to use current axes.
        grid_spacing : `~astropy.units.Quantity`
            Spacing for longitude and latitude grid, if length two it specifies
            (lon, lat) spacing.
        annotate : `bool`
            Passing `False` disables the axes labels and the ticks on the top and right axes.
        system : str
            Coordinate system for the grid. Must be 'stonyhurst' or 'carrington'.
        kwargs :
            Additional keyword arguments are passed to `~sunpy.visualization.wcsaxes_compat.wcsaxes_heliographic_overlay`.

        Returns
        -------
        overlay: `~astropy.visualization.wcsaxes.CoordinatesMap`
            The wcsaxes coordinate overlay instance.

        Notes
        -----
        Keyword arguments are passed onto the `sunpy.visualization.wcsaxes_compat.wcsaxes_heliographic_overlay` function.
        """
        axes = self._check_axes(axes)
        return wcsaxes_compat.wcsaxes_heliographic_overlay(axes,
                                                           grid_spacing=grid_spacing,
                                                           annotate=annotate,
                                                           obstime=self.reference_date,
                                                           rsun=self.rsun_meters,
                                                           observer=self.observer_coordinate,
                                                           system=system,
                                                           **kwargs)

    def draw_limb(self, axes=None, *, resolution=1000, **kwargs):
        """
        Draws the solar limb as seen by the map's observer.

        The limb is a circle for only the simplest plots. If the coordinate frame of
        the limb is different from the coordinate frame of the plot axes, not only
        may the limb not be a true circle, a portion of the limb may be hidden from
        the observer. In that case, the circle is divided into visible and hidden
        segments, represented by solid and dotted lines, respectively.

        Parameters
        ----------
        axes : `~matplotlib.axes` or ``None``
            Axes to plot limb on or ``None`` to use current axes.
        resolution : `int`
            The number of points to use to represent the limb.

        Returns
        -------
        visible : `~matplotlib.patches.Polygon` or `~matplotlib.patches.Circle`
            The patch added to the axes for the visible part of the limb (i.e., the
            "near" side of the Sun).
        hidden : `~matplotlib.patches.Polygon` or None
            The patch added to the axes for the hidden part of the limb (i.e., the
            "far" side of the Sun).

        Notes
        -----
        Keyword arguments are passed onto the patches.

        If the limb is a true circle, ``visible`` will instead be
        `~matplotlib.patches.Circle` and ``hidden`` will be ``None``. If there
        are no visible points (e.g., on a synoptic map any limb is fully)
        visible ``hidden`` will be ``None``.

        To avoid triggering Matplotlib auto-scaling, these patches are added as
        artists instead of patches. One consequence is that the plot legend is not
        populated automatically when the limb is specified with a text label. See
        :ref:`sphx_glr_gallery_text_labels_and_annotations_custom_legends.py` in
        the Matplotlib documentation for examples of creating a custom legend.
        """
        # Put imports here to reduce sunpy.map import time
        import sunpy.visualization.drawing

        axes = self._check_axes(axes)
        return sunpy.visualization.drawing.limb(
            axes,
            self.observer_coordinate,
            resolution=resolution,
            rsun=self.rsun_meters,
            **kwargs
        )

    def draw_extent(self, *, axes=None, **kwargs):
        """
        Draw the extent of the map onto a given axes.

        Parameters
        ----------
        axes : `matplotlib.axes.Axes`, optional
            The axes to plot the extent on, or "None" to use current axes.

        Returns
        -------
        visible : `~matplotlib.patches.Polygon`
            The patch added to the axes for the visible part of the WCS extent.
        hidden : `~matplotlib.patches.Polygon`
            The patch added to the axes for the hidden part of the WCS extent.
        """
        # Put imports here to reduce sunpy.map import time
        import sunpy.visualization.drawing

        axes = self._check_axes(axes)
        return sunpy.visualization.drawing.extent(
            axes,
            self.wcs,
            **kwargs
        )

    @u.quantity_input
    def draw_quadrangle(self, bottom_left, *, width: (u.deg, u.pix) = None, height: (u.deg, u.pix) = None,
                        axes=None, top_right=None, **kwargs):
        """
        Draw a quadrangle defined in world coordinates on the plot using Astropy's
        `~astropy.visualization.wcsaxes.Quadrangle`.

        This draws a quadrangle that has corners at ``(bottom_left, top_right)``,
        and has sides aligned with the coordinate axes of the frame of ``bottom_left``,
        which may be different from the coordinate axes of the map.

        If ``width`` and ``height`` are specified, they are respectively added to the
        longitude and latitude of the ``bottom_left`` coordinate to calculate a
        ``top_right`` coordinate.

        Parameters
        ----------
        bottom_left : `~astropy.coordinates.SkyCoord` or `~astropy.units.Quantity`
            The bottom-left coordinate of the rectangle. If a `~astropy.coordinates.SkyCoord` it can
            have shape ``(2,)`` and simultaneously define ``top_right``. If specifying
            pixel coordinates it must be given as an `~astropy.units.Quantity`
            object with pixel units (e.g., ``pix``).
        top_right : `~astropy.coordinates.SkyCoord` or `~astropy.units.Quantity`, optional
            The top-right coordinate of the quadrangle. If ``top_right`` is
            specified ``width`` and ``height`` must be omitted.
        width : `astropy.units.Quantity`, optional
            The width of the quadrangle. Required if ``top_right`` is omitted.
        height : `astropy.units.Quantity`
            The height of the quadrangle. Required if ``top_right`` is omitted.
        axes : `matplotlib.axes.Axes`
            The axes on which to plot the quadrangle. Defaults to the current
            axes.

        Returns
        -------
        quad : `~astropy.visualization.wcsaxes.Quadrangle`
            The added patch.

        Notes
        -----
        Extra keyword arguments to this function are passed through to the
        `~astropy.visualization.wcsaxes.Quadrangle` instance.

        Examples
        --------
        .. minigallery:: sunpy.map.GenericMap.draw_quadrangle
        """
        axes = self._check_axes(axes)

        if isinstance(bottom_left, u.Quantity):
            anchor, _, top_right, _ = self._parse_submap_quantity_input(bottom_left, top_right, width, height)
            width, height = top_right - anchor
            transform = axes.get_transform(self.wcs if self.wcs is not axes.wcs else 'pixel')
            kwargs.update({"vertex_unit": u.pix})

        else:
            bottom_left, top_right = get_rectangle_coordinates(
                bottom_left, top_right=top_right, width=width, height=height)

            width = Longitude(top_right.spherical.lon - bottom_left.spherical.lon)
            height = top_right.spherical.lat - bottom_left.spherical.lat
            anchor = self._get_lon_lat(bottom_left)
            transform = axes.get_transform(bottom_left.replicate_without_data())

        kwergs = {
            "transform": transform,
            "edgecolor": "white",
            "fill": False,
        }
        kwergs.update(kwargs)
        quad = Quadrangle(anchor, width, height, **kwergs)
        axes.add_patch(quad)
        return quad

    def _process_levels_arg(self, levels):
        """
        Accept a percentage or dimensionless or map unit input for contours.
        """
        levels = np.atleast_1d(levels)
        if not hasattr(levels, 'unit'):
            if self.unit is None:
                # No map units, so allow non-quantity through
                return levels
            else:
                raise TypeError("The levels argument has no unit attribute, "
                                "it should be an Astropy Quantity object.")

        if levels.unit == u.percent:
            return 0.01 * levels.to_value('percent') * np.nanmax(self.data)
        elif self.unit is not None:
            return levels.to_value(self.unit)
        elif levels.unit.is_equivalent(u.dimensionless_unscaled):
            # Handle case where map data has no units
            return levels.to_value(u.dimensionless_unscaled)
        else:
            # Map data has no units, but levels doesn't have dimensionless units
            raise u.UnitsError("This map has no unit, so levels can only be specified in percent "
                               "or in u.dimensionless_unscaled units.")


    def _update_contour_args(self, contour_args):
        """
        Updates ``contour_args`` with values from ``plot_settings``.

        Parameters
        ----------
        contour_args : dict
            A dictionary of arguments to be used for contour plotting.

        Returns
        -------
        dict
            The updated ``contour_args`` dictionary.

        Notes
        -----
        - 'cmap': Set to `None` to avoid the error "ValueError: Either colors or cmap must be None".
        - 'interpolation': Removed because Matplotlib's contour function raises the warning
        "The following kwargs were not used by contour: 'interpolation'".
        - 'origin': If `'origin': 'lower'` is present, it is replaced with `'origin': None`,
        as `None` is the default value for Matplotlib's contour plots.
        """
        plot_settings = self.plot_settings.copy()
        contour_args_copy = contour_args.copy()
        contour_args.update(plot_settings)
        # Define default settings for normal plots and contour-specific updates
        original_plot_defaults = {
            'origin': 'lower',
        }
        default_contour_param = {
            'origin': None,
        }
        # Replace conflicting settings with contour defaults
        for key in original_plot_defaults:
            if key in contour_args and contour_args[key] == original_plot_defaults[key]:
                contour_args[key] = default_contour_param[key]
        # 'cmap' cannot be used for contour plots when levels are not None,
        # which is the case in composite maps.
        contour_args['cmap'] = None
        # custom 'norm' cannot be passed through plot_settings
        contour_args['norm'] = None
        # If 'draw_contour' is used, setting 'norm' and 'cmap' to None ensures the method arguments are applied.
        contour_args.update(contour_args_copy)
        contour_args.pop('interpolation')
        return contour_args


    def draw_contours(self, levels, axes=None, *, fill=False, **contour_args):
        """
        Draw contours of the data.

        Parameters
        ----------
        levels : `~astropy.units.Quantity`
            A list of numbers indicating the contours to draw. These are given
            as a percentage of the maximum value of the map data, or in units
            equivalent to the `~sunpy.map.GenericMap.unit` attribute.
        axes : `matplotlib.axes.Axes`
            The axes on which to plot the contours. Defaults to the current
            axes.
        fill : `bool`, optional
            Determines the style of the contours:
            - If `False` (default), contours are drawn as lines using :meth:`~matplotlib.axes.Axes.contour`.
            - If `True`, contours are drawn as filled regions using :meth:`~matplotlib.axes.Axes.contourf`.

        Returns
        -------
        cs : `list`
            The `~matplotlib.contour.QuadContourSet` object, after it has been added to
            ``axes``.

        Notes
        -----
        Extra keyword arguments to this function are passed through to the
        corresponding matplotlib method.
        """
        contour_args = self._update_contour_args(contour_args)

        axes = self._check_axes(axes)
        levels = self._process_levels_arg(levels)

        # Pixel indices
        y, x = np.indices(self.data.shape)

        # Prepare a local variable in case we need to mask values
        data = self.data

        # Transform the indices if plotting to a different WCS
        # We do this instead of using the `transform` keyword argument so that Matplotlib does not
        # get confused about the bounds of the contours
        if self.wcs is not axes.wcs:
            if "transform" in contour_args:
                transform_orig = contour_args.pop("transform")
            else:
                transform_orig = axes.get_transform(self.wcs)
            transform = transform_orig - axes.transData  # pixel->pixel transform
            x_1d, y_1d = transform.transform(np.stack([x.ravel(), y.ravel()]).T).T
            x, y = np.reshape(x_1d, x.shape), np.reshape(y_1d, y.shape)

            # Mask out the data array anywhere the coordinate arrays are not finite
            data = np.ma.array(data, mask=~np.logical_and(np.isfinite(x), np.isfinite(y)))

        if fill:
            # Ensure we have more than one level if fill is True
            if len(levels) == 1:
                max_val = np.nanmax(self.data)
                # Ensure the existing level is less than max_val
                if levels[0] < max_val:
                    levels = np.append(levels, max_val)
                else:
                    raise ValueError(
                        f"The provided level ({levels[0]}) is not smaller than the maximum data value ({max_val}). "
                        "Contour level must be smaller than the maximum data value to use `fill=True`.")
            cs = axes.contourf(x, y, data, levels, **contour_args)
        else:
            cs = axes.contour(x, y, data, levels, **contour_args)
        return cs

    @peek_show
    def peek(self, draw_limb=False, draw_grid=False,
             colorbar=True, **matplot_args):
        """
        Displays a graphical overview of the data in this object for user evaluation.
        For the creation of plots, users should instead use the `~sunpy.map.GenericMap.plot`
        method and Matplotlib's pyplot framework.

        Parameters
        ----------
        draw_limb : bool
            Whether the solar limb should be plotted.
        draw_grid : bool or `~astropy.units.Quantity`
            Whether solar meridians and parallels are plotted.
            If `~astropy.units.Quantity` then sets degree difference between
            parallels and meridians.
        colorbar : bool
            Whether to display a colorbar next to the plot.
        **matplot_args : dict
            Matplotlib Any additional imshow arguments that should be used
            when plotting.
        """
        figure = plt.figure()
        axes = wcsaxes_compat.gca_wcs(self.wcs)

        im = self.plot(axes=axes, **matplot_args)

        grid_spacing = None
        # Handle case where draw_grid is actually the grid sapcing
        if isinstance(draw_grid, u.Quantity):
            grid_spacing = draw_grid
            draw_grid = True
        elif not isinstance(draw_grid, bool):
            raise TypeError("draw_grid should be a bool or an astropy Quantity.")

        if colorbar:
            if draw_grid:
                pad = 0.12  # Pad to compensate for ticks and axes labels
            else:
                pad = 0.05  # Default value for vertical colorbar
            colorbar_label = str(self.unit) if self.unit is not None else ""
            figure.colorbar(im, pad=pad).set_label(colorbar_label,
                                                   rotation=0, labelpad=-50, y=-0.02, size=12)

        if draw_limb:
            self.draw_limb(axes=axes)

        if draw_grid:
            if grid_spacing is None:
                self.draw_grid(axes=axes)
            else:
                self.draw_grid(axes=axes, grid_spacing=grid_spacing)

        return figure

    @u.quantity_input
    def plot(self, *, annotate=True, axes=None, title=True, autoalign=True,
             clip_interval: u.percent = None, **imshow_kwargs):
        """
        Plots the map using matplotlib.

        By default, the map's pixels will be drawn in an coordinate-aware fashion, even
        when the plot axes are a different projection or a different coordinate frame.
        See the ``autoalign`` keyword and the notes below.

        Parameters
        ----------
        annotate : `bool`, optional
            If `True`, the data is plotted at its natural scale; with
            title and axis labels.
        axes : `~matplotlib.axes.Axes` or None
            If provided the image will be plotted on the given axes. Else the
            current Matplotlib axes will be used.
        title : `str`, `bool`, optional
            The plot title. If `True`, uses the default title for this map.
        clip_interval : two-element `~astropy.units.Quantity`, optional
            If provided, the data will be clipped to the percentile interval bounded by the two
            numbers.
        autoalign : `bool` or `str`, optional
            If other than `False`, the plotting accounts for any difference between the
            WCS of the map and the WCS of the `~astropy.visualization.wcsaxes.WCSAxes`
            axes (e.g., a difference in rotation angle). The options are:

            * ``"mesh"``, which draws a mesh of the individual map pixels
            * ``"image"``, which draws the map as a single (warped) image
            * `True`, which automatically determines whether to use ``"mesh"`` or ``"image"``

        **imshow_kwargs : `dict`
            Any additional arguments are passed to :meth:`~matplotlib.axes.Axes.imshow`
            or :meth:`~matplotlib.axes.Axes.pcolormesh`.

        Examples
        --------
        >>> # Simple Plot with color bar
        >>> aia.plot()   # doctest: +SKIP
        >>> plt.colorbar()   # doctest: +SKIP
        >>> # Add a limb line and grid
        >>> aia.plot()   # doctest: +SKIP
        >>> aia.draw_limb()   # doctest: +SKIP
        >>> aia.draw_grid()   # doctest: +SKIP

        Notes
        -----
        The ``autoalign`` functionality can be intensive to render. If the plot is to
        be interactive, the alternative approach of preprocessing the map to match the
        intended axes (e.g., through rotation or reprojection) will result in better
        plotting performance.

        The ``autoalign='image'`` approach is usually faster than the
        ``autoalign='mesh'`` approach, but is not as reliable, depending on the
        specifics of the map.  If parts of the map cannot be plotted, a warning is
        emitted.  If the entire map cannot be plotted, an error is raised.

        When combining ``autoalign`` functionality with
        `~sunpy.coordinates.Helioprojective` coordinates, portions of the map that are
        beyond the solar disk may not appear.  To preserve the off-disk parts of the
        map, using the `~sunpy.coordinates.SphericalScreen` context manager may be
        appropriate.
        """
        if autoalign == 'pcolormesh':
            warn_deprecated("Specifying `autoalign='pcolormesh'` is deprecated as of 7.0. "
                            "Specify `autoalign='mesh'` instead.")
            autoalign = 'mesh'

        # Set the default approach to autoalignment
        allowed_autoalign = [False, True, 'mesh', 'image']
        if autoalign not in allowed_autoalign:
            raise ValueError(f"The value for `autoalign` must be one of {allowed_autoalign}.")

        axes = self._check_axes(axes, warn_different_wcs=autoalign is False)

        # Normal plot
        plot_settings = copy.deepcopy(self.plot_settings)
        if 'title' in plot_settings:
            plot_settings_title = plot_settings.pop('title')
        else:
            plot_settings_title = self.latex_name

        # Anything left in plot_settings is given to imshow
        imshow_args = plot_settings
        if annotate:
            if title is True:
                title = plot_settings_title

            if title:
                axes.set_title(title)

            # WCSAxes has unit identifiers on the tick labels, so no need
            # to add unit information to the label
            ctype = axes.wcs.wcs.ctype
            axes.coords[0].set_axislabel(axis_labels_from_ctype(ctype[0], None))
            axes.coords[1].set_axislabel(axis_labels_from_ctype(ctype[1], None))

        # Take a deep copy here so that a norm in imshow_kwargs doesn't get modified
        # by setting it's vmin and vmax
        imshow_args.update(copy.deepcopy(imshow_kwargs))

        if clip_interval is not None:
            vmin, vmax = _clip_interval(self.data, clip_interval)
            imshow_args['vmin'] = vmin
            imshow_args['vmax'] = vmax

        if (norm := imshow_args.get('norm', None)) is not None:
            _handle_norm(norm, imshow_args)

        if self.mask is None:
            data = self.data
        else:
            data = np.ma.array(np.asarray(self.data), mask=self.mask)

        # Disable autoalignment if it is not necessary
        # TODO: revisit tolerance value
        if autoalign is True and axes.wcs.wcs.compare(self.wcs.wcs, tolerance=0.01):
            autoalign = False

        if autoalign in {True, 'image'}:
            ny, nx = self.data.shape
            pixel_perimeter = grid_perimeter(nx, ny) - 0.5

            transform = axes.get_transform(self.wcs) - axes.transData
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=SunpyUserWarning)
                data_perimeter = transform.transform(pixel_perimeter)
            if not np.all(np.isfinite(data_perimeter)):
                if autoalign == 'image':
                    raise RuntimeError("Cannot draw an autoaligned image at all due to its coordinates. "
                                       "Try specifying autoalign=mesh.")
                autoalign = 'mesh'
            else:
                min_x, min_y = np.min(data_perimeter, axis=0)
                max_x, max_y = np.max(data_perimeter, axis=0)

                data_corners = data_perimeter[[0, ny, nx + ny, nx + 2*ny], :]
                if not (np.allclose([min_x, min_y], np.min(data_corners, axis=0))
                        and np.allclose([max_x, max_y], np.max(data_corners, axis=0))):
                    if autoalign == 'image':
                        warn_user("Cannot draw all of the autoaligned image due to the warping required. "
                                  "Specifying autoalign=mesh is recommended.")
                    else:
                        autoalign = 'mesh'
            if autoalign == 'mesh':
                log.info("Using mesh-based autoalignment")
            elif autoalign is True:
                log.info("Using image-based autoalignment")
                autoalign = 'image'

        if autoalign == 'image':
            # Draw the image, but revert to the prior data limits because matplotlib does not account for the transform
            old_datalim = copy.deepcopy(axes.dataLim)
            ret = axes.imshow(data, transform=transform + axes.transData, **imshow_args)
            axes.dataLim = old_datalim

            # Update the data limits based on the transformed perimeter
            ret.sticky_edges.x[:] = [min_x, max_x]
            ret.sticky_edges.y[:] = [min_y, max_y]
            axes.update_datalim([(min_x, min_y), (max_x, max_y)])
            axes.autoscale(enable=None)

            # Clip the drawn image based on the transformed perimeter
            path = matplotlib.path.Path(data_perimeter)
            ret.set_clip_path(path, axes.transData)
        elif autoalign == 'mesh':
            # We have to handle an `aspect` keyword separately
            axes.set_aspect(imshow_args.get('aspect', 1))

            # pcolormesh does not do interpolation
            if imshow_args.get('interpolation', None) not in [None, 'none', 'nearest']:
                warn_user("The interpolation keyword argument is ignored when using autoalign "
                          "functionality.")

            # Set the zorder to be 0 so that it is treated like an image in ordering
            imshow_args.setdefault('zorder', 0)

            # Remove imshow keyword arguments that are not accepted by pcolormesh
            for item in ['aspect', 'extent', 'interpolation', 'origin']:
                if item in imshow_args:
                    del imshow_args[item]

            # The quadrilaterals of pcolormesh can slightly overlap, which creates the appearance
            # of a grid pattern when alpha is not 1. These settings minimize the overlap.
            if imshow_args.get('alpha', 1) != 1:
                imshow_args.setdefault('antialiased', True)
                imshow_args.setdefault('linewidth', 0)

            # Create a lookup table for the transformed data corners for matplotlib to use
            transform = _PrecomputedPixelCornersTransform(axes, self.wcs)

            # Define the mesh in data coordinates in case the transformation results in NaNs
            ret = axes.pcolormesh(transform.data_x, transform.data_y, data,
                                  shading='flat',
                                  transform=transform + axes.transData,
                                  **imshow_args)

            # Calculate the bounds of the mesh in the pixel space of the axes
            good = np.logical_and(np.isfinite(transform.axes_x), np.isfinite(transform.axes_y))
            good_x, good_y = transform.axes_x[good], transform.axes_y[good]
            min_x, max_x = np.min(good_x), np.max(good_x)
            min_y, max_y = np.min(good_y), np.max(good_y)

            # Update the plot limits
            ret.sticky_edges.x[:] = [min_x, max_x]
            ret.sticky_edges.y[:] = [min_y, max_y]
            axes.update_datalim([(min_x, min_y), (max_x, max_y)])
        else:
            ret = axes.imshow(data, **imshow_args)

        wcsaxes_compat.default_wcs_grid(axes)

        # Set current axes/image if pyplot is being used (makes colorbar work)
        for i in plt.get_fignums():
            if axes in plt.figure(i).axes:
                plt.sca(axes)
                plt.sci(ret)

        return ret

    @deprecated(since="6.1", alternative="sunpy.map.GenericMap.find_contours")
    def contour(self, level, **kwargs):
        """
        Returns coordinates of the contours for a given level value.

        For details of the contouring algorithm see `skimage.measure.find_contours`.

        Parameters
        ----------
        level : float, `~astropy.units.Quantity`
            Value along which to find contours in the array. If the map unit attribute
            is not `None`, this must be a `~astropy.units.Quantity` with units
            equivalent to the map data units.
        kwargs :
            Additional keyword arguments are passed to `skimage.measure.find_contours`.

        Returns
        -------
        contours: list of (n,2) `~astropy.coordinates.SkyCoord`
            Coordinates of each contour.

        Examples
        --------
        >>> import astropy.units as u
        >>> import sunpy.map
        >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
        >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
        >>> contours = aia.contour(50000 * u.DN)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
        >>> contours[0]  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
        <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.880, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.880, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
        (-0.00406429, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
        [(719.59798458, -352.60839064), (717.19243987, -353.75348121),
         (715.8820808 , -354.75140718), (714.78652558, -355.05102034),
         (712.68209174, -357.14645009), (712.68639008, -359.54923801),
         (713.14112796, -361.95311455), (714.76598031, -363.53013567),
         (717.17229147, -362.06880784), (717.27714042, -361.9631112 ),
         (718.43620686, -359.56313541), (718.8672722 , -357.1614    ),
         (719.58811599, -356.68119768), (721.29217122, -354.76448374),
         (719.59798458, -352.60839064)]>

        See Also
        --------
        skimage.measure.find_contours
        """
        from skimage import measure

        level = self._process_levels_arg(level)
        if level.size != 1:
            raise ValueError("level must be a single scalar value")
        else:
            # _process_levels_arg converts level to a 1D array, but
            # find_contours expects a scalar below
            level = level[0]

        contours = measure.find_contours(self.data, level=level, **kwargs)
        contours = [self.wcs.array_index_to_world(c[:, 0], c[:, 1]) for c in contours]
        return contours

    def find_contours(self, level, method='contourpy', **kwargs):
        """
        Returns coordinates of the contours for a given level value.

        For details of the contouring algorithm, see :func:`contourpy.contour_generator` or :func:`skimage.measure.find_contours`.

        Parameters
        ----------
        level : float, astropy.units.Quantity
            Value along which to find contours in the array. If the map unit attribute
            is not `None`, this must be a `~astropy.units.Quantity` with units
            equivalent to the map data units.
        method : {'contourpy', 'skimage'}
            Determines which contouring method is used and should
            be specified as either 'contourpy' or 'skimage'.
            Defaults to 'contourpy'.
        kwargs :
            Additional keyword arguments passed to either :func:`contourpy.contour_generator`
            or :func:`skimage.measure.find_contours`, depending on the value of the ``method`` argument.

        Returns
        -------
        contours: list of (n,2) `~astropy.coordinates.SkyCoord`
            Coordinates of each contour.

        Examples
        --------
        >>> import astropy.units as u
        >>> import sunpy.map
        >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
        >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA +IGNORE_WARNINGS
        >>> contours = aia.find_contours(50000 * u.DN, method='contourpy')  # doctest: +REMOTE_DATA
        >>> contours[0]  # doctest: +REMOTE_DATA
        <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.880, rsun=696000.0 km, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.880, rsun=696000.0 km): (lon, lat, radius) in (deg, deg, m)
            (-0.00406429, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
            [(713.14112796, -361.95311455), (714.76598031, -363.53013567),
             (717.17229147, -362.06880784), (717.27714042, -361.9631112 ),
             (718.43620686, -359.56313541), (718.8672722 , -357.1614    ),
             (719.58811599, -356.68119768), (721.29217122, -354.76448374),
             (722.00110323, -352.46446792), (722.08933899, -352.36363319),
             (722.00223989, -351.99536019), (721.67724425, -352.36263712),
             (719.59798458, -352.60839064), (717.19243987, -353.75348121),
             (715.8820808 , -354.75140718), (714.78652558, -355.05102034),
             (712.68209174, -357.14645009), (712.68639008, -359.54923801),
             (713.14112796, -361.95311455)]>

        See Also
        --------
        :func:`contourpy.contour_generator`
        :func:`skimage.measure.find_contours`
        """
        level = self._process_levels_arg(level)
        if level.size != 1:
            raise ValueError("level must be a single scalar value")
        else:
            # _process_levels_arg converts level to a 1D array, but
            # find_contours expects a scalar below
            level = level[0]

        if method == 'contourpy':
            from contourpy import contour_generator

            gen = contour_generator(z=self.data, **kwargs)
            contours = gen.lines(level)
            contours = [self.wcs.array_index_to_world(c[:, 1], c[:, 0]) for c in contours]
        elif method == 'skimage':
            from skimage import measure

            contours = measure.find_contours(self.data, level=level, **kwargs)
            contours = [self.wcs.array_index_to_world(c[:, 0], c[:, 1]) for c in contours]
        else:
            raise ValueError(f"Unknown method '{method}'. Use 'contourpy' or 'skimage'.")

        return contours

    def _check_axes(self, axes, warn_different_wcs=False):
        """
        - If axes is None, get the current Axes object.
        - Error if not a WCSAxes.
        - Return axes.

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            Axes to validate.
        warn_different_wcs : bool
            If `True`, warn if the Axes WCS is different from the Map WCS. This is only used for
            `.plot()`, and can be removed once support is added for plotting a map on a different
            WCSAxes.
        """
        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        if not wcsaxes_compat.is_wcsaxes(axes):
            raise TypeError("The axes need to be an instance of WCSAxes. "
                            "To fix this pass set the `projection` keyword "
                            "to this map when creating the axes.")
        elif warn_different_wcs and not axes.wcs.wcs.compare(self.wcs.wcs, tolerance=0.01):
            warn_user('The map world coordinate system (WCS) is different from the axes WCS. '
                      'The map data axes may not correctly align with the coordinate axes. '
                      'To automatically transform the data to the coordinate axes, specify '
                      '`autoalign=True`.')

        return axes

    def reproject_to(self, target_wcs, *, algorithm='interpolation', return_footprint=False,
                     auto_extent: Literal[None, 'corners', 'edges', 'all'] = None,
                     **reproject_args):
        """
        Reproject the map to a different world coordinate system (WCS)

        Additional keyword arguments are passed through to the reprojection function.

        This method **does not** preserve dask arrays.

        Parameters
        ----------
        target_wcs : `dict` or `~astropy.wcs.WCS`
            The destination FITS WCS header or WCS instance
        algorithm : `str`
            One of the supported `reproject` algorithms (see below)
        return_footprint : `bool`
            If ``True``, the footprint is returned in addition to the new map.
            Defaults to ``False``.
        auto_extent : ``"all"``, ``"edges"``, ``"corners"``, or ``None``
            If ``None``, the extent of the reprojected map comes from the target WCS.
            If not ``None``, the extent of the reprojected map is automatically
            determined by ensuring that all of the pixels, just the edges, or just the
            corners of this map are in the contained within the extent.  Compared to the
            target WCS, the extent will be shifted/expanded/cropped by an integer number
            of pixels.
            Defaults to ``None``.

        Returns
        -------
        outmap : `~sunpy.map.GenericMap`
            The reprojected map
        footprint : `~numpy.ndarray`
            Footprint of the input arary in the output array. Values of 0 indicate no
            coverage or valid values in the input image, while values of 1 indicate
            valid values. Intermediate values indicate partial coverage.
            Only returned if ``return_footprint`` is ``True``.

        Notes
        -----
        The reprojected map does not preserve any metadata beyond the WCS-associated
        metadata.

        The supported `reproject` algorithms are:

        * 'interpolation' for :func:`~reproject.reproject_interp`
        * 'adaptive' for :func:`~reproject.reproject_adaptive`
        * 'exact' for :func:`~reproject.reproject_exact`

        See the respective documentation for these functions for additional keyword
        arguments that are allowed.

        Of the options for the automatic determination of the reprojected extent, both
        ``"edges"`` and ``"corners"`` will perform the calculation faster than
        ``"all"``, but at the risk of potentially not including the entire reprojected
        map.

        .. minigallery:: sunpy.map.GenericMap.reproject_to
        """
        if not isinstance(target_wcs, astropy.wcs.WCS):
            target_wcs = astropy.wcs.WCS(target_wcs)

        # Select the desired reprojection algorithm
        functions = {'interpolation': reproject.reproject_interp,
                     'adaptive': reproject.reproject_adaptive,
                     'exact': reproject.reproject_exact}
        if algorithm not in functions:
            raise ValueError(f"The specified algorithm must be one of: {list(functions.keys())}")
        func = functions[algorithm]

        if auto_extent not in ['all', 'edges', 'corners', None]:
            raise ValueError("The allowed options for `auto_extent` are 'all', 'edges', 'corners', or None.")
        if auto_extent is not None:
            left, right, bottom, top = extent_in_other_wcs(self.wcs, target_wcs, original_shape=self.data.shape,
                                                            method=auto_extent, integers=True)
            target_wcs.wcs.crpix -= [left, bottom]
            target_wcs.pixel_shape = [right - left + 1, top - bottom + 1]

        # reproject does not automatically grab the array shape from the WCS instance
        if target_wcs.array_shape is not None:
            reproject_args.setdefault('shape_out', target_wcs.array_shape)

        # Reproject the array
        output_array = func(self, target_wcs, return_footprint=return_footprint, **reproject_args)
        if return_footprint:
            output_array, footprint = output_array

        # Create and return a new GenericMap
        outmap = GenericMap(output_array, target_wcs.to_header(),
                            plot_settings=self.plot_settings)

        # Check rsun mismatch
        if self.rsun_meters != outmap.rsun_meters:
            warn_user("rsun mismatch detected: "
                      f"{self.name}.rsun_meters={self.rsun_meters}; {outmap.name}.rsun_meters={outmap.rsun_meters}. "
                      "This might cause unexpected results during reprojection.")

        if return_footprint:
            return outmap, footprint
        return outmap


GenericMap.__doc__ = fix_duplicate_notes(_notes_doc, GenericMap.__doc__)


class InvalidHeaderInformation(ValueError):
    """Exception to raise when an invalid header tag value is encountered for a
    FITS/JPEG 2000 file."""
