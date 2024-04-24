"""
Map is a generic Map class from which all other Map classes inherit from.
"""
import copy
import html
import inspect
import numbers
import textwrap
import itertools
import webbrowser
from numbers import Integral
from tempfile import NamedTemporaryFile

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backend_bases import FigureCanvasBase
from matplotlib.figure import Figure

try:
    from dask.array import Array as DaskArray
    DASK_INSTALLED = True
except ImportError:
    DASK_INSTALLED = False

import astropy.units as u
import astropy.wcs
from astropy.coordinates import BaseCoordinateFrame, SkyCoord, UnitSphericalRepresentation
from astropy.utils.metadata import MetaData
from astropy.visualization import HistEqStretch, ImageNormalize
from astropy.visualization.wcsaxes import WCSAxes

import sunpy.coordinates.wcs_utils
from ndcube import NDCube
from ndcube.wcs.tools import unwrap_wcs_to_fitswcs
from sunpy import config, log
# The next two are not used but are called to register functions with external modules
from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.image.resample import resample as sunpy_image_resample
from sunpy.image.resample import reshape_image_to_4d_superpixel
from sunpy.image.transform import _get_transform_method, _rotation_function_names, affine_transform
from sunpy.io._file_tools import write_file
from sunpy.map.mixins.mapdeprecate import MapDeprecateMixin
from sunpy.map.mixins.mapmeta import MapMetaMixin
from sunpy.util import MetaDict
from sunpy.util.decorators import ACTIVE_CONTEXTS, add_common_docstring, deprecated
from sunpy.util.exceptions import warn_user
from sunpy.util.functools import seconddispatch
from sunpy.util.util import _figure_to_base64, fix_duplicate_notes
from sunpy.visualization.plotter.mpl_plotter import MapPlotter

TIME_FORMAT = config.get("general", "time_format")
_NUMPY_COPY_IF_NEEDED = False if np.__version__.startswith("1.") else None

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

__all__ = ['GenericMap']


class GenericMap(MapDeprecateMixin, MapMetaMixin, NDCube):
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
    Reference Date:              2011-06-07 06:33:02
    Exposure Time:               0.234256 s
    Pixel Dimensions:            [1024. 1024.]
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

    def __init__(self, data, *, wcs=None, uncertainty=None, mask=None, meta,
                 unit=None, copy=False, plot_settings=None, **kwargs):
        # Setup some attributes
        self._metadata_validated = False
        self._nickname = None
        # These are placeholders for default attributes, which are only set
        # once if their data isn't present in the map metadata.
        self._default_time = None
        self._default_dsun = None
        self._default_carrington_longitude = None
        self._default_heliographic_latitude = None
        self._default_heliographic_longitude = None

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

        params = list(inspect.signature(NDCube).parameters)
        ndcube_kwargs = {x: kwargs.pop(x) for x in params & kwargs.keys()}
        super().__init__(data, wcs=None, uncertainty=uncertainty, mask=mask,
                         meta=MetaDict(meta), unit=unit, copy=copy,
                         **ndcube_kwargs)
        # NDData.__init__ sets self.wcs before it sets self.meta as our wcs
        # setter needs self.meta to exist we call the parent __init__ with
        # wcs=None and then set self.wcs so that meta is already set before the
        # wcs setter is run with the "real" wcs.
        self.wcs = wcs

        # Validate header
        # TODO: This should be a function of the header, not of the map
        self._validate_meta()

        self.plotter = MapPlotter
        if plot_settings:
            self.plotter.plot_settings.update(plot_settings)

    @property
    # @deprecated('6.0', alternative='plotter.plot_settings')
    def plot_settings(self):
        return self.plotter.plot_settings

    @plot_settings.setter
    # @deprecated('6.0', alternative='plotter.plot_settings')
    def plot_settings(self, value):
        self.plotter.plot_settings = value

    def __getitem__(self, key):
        def format_slice(key):
            if not isinstance(key, slice):
                return f"{key}"
            return f"{key.start or ''}:{key.stop or ''}{f':{key.step}' if key.step else ''}"

        if not isinstance(key, tuple):
            key = (key,)

        if any(intidx := list(map(lambda x: isinstance(x, Integral), key))):
            strslice = ", ".join([f"{k}:{k+1}" if isint else format_slice(k)
                                  for k, isint in zip(key, intidx)])
            raise TypeError(
                "It is not possible to slice a map with an integer as it will "
                "reduce the number of data dimensions by one.\nIn order to "
                "apply the same slice without dropping a dimension do "
                f"mymap[{strslice}]."
            )
        return super().__getitem__(key)


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
                   Reference Date:\t\t {ref_date}
                   Exposure Time:\t\t {dt}
                   Pixel Dimensions:\t\t {dim}
                   Coordinate System:\t {coord}
                   Scale:\t\t\t {scale}
                   Reference Pixel:\t {refpix}
                   Reference Coord:\t {refcoord}\
                   """).format(obs=self.observatory, inst=self.instrument, det=self.detector,
                               meas=measurement, wave=wave,
                               date=self.date.strftime(TIME_FORMAT),
                               ref_date=self.reference_date.strftime(TIME_FORMAT),
                               dt=dt,
                               dim=u.Quantity(self.shape[::-1]),
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

    def _new_instance(self, data=None, meta=None, plot_settings=None, **kwargs):
        """
        Instantiate a new instance of this class using given data.
        This is a shortcut for ``type(self)(data, meta, plot_settings)``.
        """
        if meta is None:
            meta = copy.deepcopy(getattr(self, 'meta'))
        if new_unit := kwargs.get('unit', None):
            meta['bunit'] = self._parse_fits_unit(new_unit).to_string()
        # NOTE: wcs=None is explicitly passed here because the wcs of a map is
        # derived from the information in the metadata.
        new_map = super()._new_instance(data=data, meta=meta, wcs=None, **kwargs)
        # plot_settings are set explicitly here as some map sources
        # explicitly set some of the plot_xsettings in the constructor
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
        new_map = self._new_instance(data=self.data, meta=new_meta, plot_settings=self.plotter.plot_settings)

        return new_map

    @property
    # TODO FIX THIS
    #@cached_property_based_on('_meta_hash')
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
        self._validate_meta()
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

    @wcs.setter
    def wcs(self, wcs):
        """
        Map uses the meta dict as the source of truth.
        When setting the wcs of the map we convert it to a header and then
        update the header of the map.
        """
        if wcs is None:
            return
        # Unwrap any wrapper classes to FITS WCS
        unwrapped, _ = unwrap_wcs_to_fitswcs(wcs)
        # Convert to a header
        new_header = unwrapped.to_header()
        # wcslib ignores NAXIS
        for n in range(1, unwrapped.naxis + 1):
            new_header[f"NAXIS{n}"] = unwrapped._naxis[n-1]
        old_wcs_header = self.wcs.to_header()
        # Reduce the new header to just the keys which differ from the current WCS
        # We do this to figure out what's been changed post wcslib doing any
        # conversion (such as arcsec -> deg)
        changed_header = dict(set(new_header.items()).difference(old_wcs_header.items()))
        # If the units are different to the original metadata then we override
        # all the keys from the new WCS header to account for unit conversion.
        if any([new_header[f"CUNIT{n}"] != self.meta[f"CUNIT{n}"] for n in range(1, 3)]):
            # TODO: Do we want to make it so we don't change units?
            changed_header = new_header
        self.meta.update(MetaDict(changed_header))

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
        as IDL's congrid routine, which apparently originally came from a
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

        scale_factor_x = float(self.shape[1] / dimensions[0].value)
        scale_factor_y = float(self.shape[0] / dimensions[1].value)

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
        new_map = self._new_instance(data=new_data, meta=new_meta, plot_settings=self.plotter.plot_settings)
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
        new_map = self._new_instance(data=new_data, meta=new_meta, plot_settings=self.plotter.plot_settings)

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
        Observatory:                 SDO
        Instrument:          AIA 3
        Detector:            AIA
        Measurement:                 171.0 Angstrom
        Wavelength:          171.0 Angstrom
        Observation Date:    2011-06-07 06:33:02
        Reference Date:              2011-06-07 06:33:02
        Exposure Time:               0.234256 s
        Pixel Dimensions:            [335. 335.]
        Coordinate System:   helioprojective
        Scale:                       [2.402792 2.402792] arcsec / pix
        Reference Pixel:     [126.5 125.5] pix
        Reference Coord:     [3.22309951 1.38578135] arcsec
        array([[ 450.4546 ,  565.81494,  585.0416 , ..., 1005.28284,  977.8161 ,
                1005.28284],
            [ 474.20004,  516.1865 ,  555.7032 , ..., 1010.1449 , 1010.1449 ,
                1121.2855 ],
            [ 548.1609 ,  620.9256 ,  620.9256 , ..., 1074.4924 , 1108.4492 ,
                1069.6414 ],
            ...,
            [ 206.00058,  212.1806 ,  232.78065, ...,  622.12177,  537.6615 ,
                574.74164],
            [ 229.32516,  236.07002,  222.5803 , ...,  586.8026 ,  591.2992 ,
                728.44464],
            [ 184.20439,  187.92569,  225.1387 , ...,  649.367  ,  686.58   ,
                673.5554 ]], dtype=float32)
        >>> aia.submap([0,0]*u.pixel, top_right=[5,5]*u.pixel)   # doctest: +REMOTE_DATA
        <sunpy.map.sources.sdo.AIAMap object at ...>
        SunPy Map
        ---------
        Observatory:                 SDO
        Instrument:          AIA 3
        Detector:            AIA
        Measurement:                 171.0 Angstrom
        Wavelength:          171.0 Angstrom
        Observation Date:    2011-06-07 06:33:02
        Reference Date:              2011-06-07 06:33:02
        Exposure Time:               0.234256 s
        Pixel Dimensions:            [6. 6.]
        Coordinate System:   helioprojective
        Scale:                       [2.402792 2.402792] arcsec / pix
        Reference Pixel:     [511.5 511.5] pix
        Reference Coord:     [3.22309951 1.38578135] arcsec
        array([[-95.92475   ,   7.076416  ,  -1.9656711 ,  -2.9485066 ,
                -0.98283553,  -6.0935802 ],
            [-96.97533   ,  -5.1167884 ,   0.        ,   0.        ,
                0.9746264 ,   3.8985057 ],
            [-93.99607   ,   1.0189276 ,  -4.0757103 ,   2.0378551 ,
                -2.0378551 ,  -7.896689  ],
            [-96.97533   ,  -8.040668  ,  -2.9238791 ,  -5.1167884 ,
                -0.9746264 ,  -8.040668  ],
            [-95.92475   ,   6.028058  ,  -4.9797    ,  -1.0483578 ,
                -3.9313421 ,  -1.0483578 ],
            [-95.103004  ,   0.        ,  -4.993475  ,   0.        ,
                -4.0855703 ,  -7.03626   ]], dtype=float32)
        >>> width = 10 * u.arcsec
        >>> height = 10 * u.arcsec
        >>> aia.submap(bl, width=width, height=height)   # doctest: +REMOTE_DATA
        <sunpy.map.sources.sdo.AIAMap object at ...>
        SunPy Map
        ---------
        Observatory:                 SDO
        Instrument:          AIA 3
        Detector:            AIA
        Measurement:                 171.0 Angstrom
        Wavelength:          171.0 Angstrom
        Observation Date:    2011-06-07 06:33:02
        Reference Date:              2011-06-07 06:33:02
        Exposure Time:               0.234256 s
        Pixel Dimensions:            [5. 5.]
        Coordinate System:   helioprojective
        Scale:                       [2.402792 2.402792] arcsec / pix
        Reference Pixel:     [125.5 125.5] pix
        Reference Coord:     [3.22309951 1.38578135] arcsec
        array([[565.81494, 585.0416 , 656.4552 , 670.18854, 678.4286 ],
            [516.1865 , 555.7032 , 634.7365 , 661.90424, 587.8105 ],
            [620.9256 , 620.9256 , 654.8825 , 596.6707 , 531.18243],
            [667.5083 , 560.52094, 651.22766, 530.28534, 495.39816],
            [570.15643, 694.5542 , 653.0883 , 699.7374 , 583.11456]],
            dtype=float32)
        >>> bottom_left_vector = SkyCoord([0, 10]  * u.deg, [0, 10] * u.deg, frame='heliographic_stonyhurst')
        >>> aia.submap(bottom_left_vector)   # doctest: +REMOTE_DATA
        <sunpy.map.sources.sdo.AIAMap object at ...>
        SunPy Map
        ---------
        Observatory:                 SDO
        Instrument:          AIA 3
        Detector:            AIA
        Measurement:                 171.0 Angstrom
        Wavelength:          171.0 Angstrom
        Observation Date:    2011-06-07 06:33:02
        Reference Date:              2011-06-07 06:33:02
        Exposure Time:               0.234256 s
        Pixel Dimensions:            [70. 69.]
        Coordinate System:   helioprojective
        Scale:                       [2.402792 2.402792] arcsec / pix
        Reference Pixel:     [1.5 0.5] pix
        Reference Coord:     [3.22309951 1.38578135] arcsec
        array([[209.89908, 213.9748 , 256.76974, ..., 560.41016, 497.23666,
                584.86444],
                [237.85315, 223.74321, 258.0102 , ..., 578.5072 , 643.00977,
                560.3659 ],
                [252.67189, 219.53459, 242.31648, ..., 623.3954 , 666.8881 ,
                625.4665 ],
                ...,
                [662.12573, 690.3013 , 702.04114, ..., 464.8968 , 561.1633 ,
                676.2135 ],
                [489.49503, 542.75616, 563.0461 , ..., 667.0321 , 748.1919 ,
                748.1919 ],
                [435.59155, 455.9701 , 496.7272 , ..., 855.8992 , 789.6689 ,
                687.7761 ]], dtype=float32)
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
            new_map = self._new_instance(data=new_data, meta=new_meta, plot_settings=self.plotter.plot_settings, mask=new_mask)
            return new_map
        # Create new map with the modification
        new_map = self._new_instance(data=new_data, meta=new_meta, plot_settings=self.plotter.plot_settings)
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
        new_map = self._new_instance(data=new_data, meta=new_meta, plot_settings=self.plotter.plot_settings, mask=new_mask)
        return new_map

    @property
    def cmap(self):
        return self.plotter.cmap

    def draw_grid(self, *args, **kwargs):
        return self.plotter.draw_grid(*args, **kwargs)

    def draw_limb(self, *args, **kwargs):
        return self.plotter.draw_limb(*args, **kwargs)

    def draw_quadrangle(self, *args, **kwargs):
        return self.plotter.draw_quadrangle(*args, **kwargs)

    def _get_cmap_name(self):
        cmap_string = f"{self.observatory}{self.detector}{self.wavelength.to_value('AA'):.0f}"
        return cmap_string.lower()

    def draw_contours(self, *args, **kwargs):
        return self.plotter.draw_contours(*args, **kwargs)

    def peek(self, *args, **kwargs):
       return self.plotter.peek(*args, **kwargs)

    def plot(self, *args, **kwargs):
        return self.plotter.plot(*args, **kwargs)

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

        For details of the contouring algorithm, see :func:`contourpy.contour_generator` or :func:`contourpy.contour_generator`.

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

    def reproject_to(self, target_wcs, *, algorithm='interpolation', return_footprint=False,
                     **reproject_args):
        """
        Reproject the map to a different world coordinate system (WCS)

        .. note::
            This method requires the optional package `reproject` to be installed.

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

        .. minigallery:: sunpy.map.GenericMap.reproject_to
        """
        # Check if both context managers are active
        if ACTIVE_CONTEXTS.get('propagate_with_solar_surface', False) and ACTIVE_CONTEXTS.get('assume_spherical_screen', False):
            warn_user("Using propagate_with_solar_surface and SphericalScreen together result in loss of off-disk data.")

        try:
            import reproject
        except ImportError as exc:
            raise ImportError("This method requires the optional package `reproject`.") from exc

        if not isinstance(target_wcs, astropy.wcs.WCS):
            target_wcs = astropy.wcs.WCS(target_wcs)

        # Select the desired reprojection algorithm
        functions = {'interpolation': reproject.reproject_interp,
                     'adaptive': reproject.reproject_adaptive,
                     'exact': reproject.reproject_exact}
        if algorithm not in functions:
            raise ValueError(f"The specified algorithm must be one of: {list(functions.keys())}")
        func = functions[algorithm]

        # reproject does not automatically grab the array shape from the WCS instance
        if target_wcs.array_shape is not None:
            reproject_args.setdefault('shape_out', target_wcs.array_shape)

        # Reproject the array
        output_array = func(self, target_wcs, return_footprint=return_footprint, **reproject_args)
        if return_footprint:
            output_array, footprint = output_array

        # Create and return a new GenericMap
        outmap = GenericMap(output_array, meta=target_wcs.to_header(),
                            plot_settings=self.plotter.plot_settings)

        # Check rsun mismatch
        if self.rsun_meters != outmap.rsun_meters:
            warn_user("rsun mismatch detected: "
                      f"{self.name}.rsun_meters={self.rsun_meters}; {outmap.name}.rsun_meters={outmap.rsun_meters}. "
                      "This might cause unexpected results during reprojection.")

        if return_footprint:
            return outmap, footprint
        return outmap


GenericMap.__doc__ += textwrap.indent(_notes_doc, "    ")


class InvalidHeaderInformation(ValueError):
    """Exception to raise when an invalid header tag value is encountered for a
    FITS/JPEG 2000 file."""
