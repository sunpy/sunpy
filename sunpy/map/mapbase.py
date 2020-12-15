"""
Map is a generic Map class from which all other Map classes inherit from.
"""
import copy
import html
import textwrap
import warnings
import webbrowser
from io import BytesIO
from base64 import b64encode
from tempfile import NamedTemporaryFile
from collections import namedtuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.backend_bases import FigureCanvasBase
from matplotlib.figure import Figure

import astropy.units as u
import astropy.wcs
from astropy.coordinates import Latitude, Longitude, SkyCoord, UnitSphericalRepresentation
from astropy.nddata import NDData
from astropy.visualization import AsymmetricPercentileInterval, HistEqStretch, ImageNormalize
from astropy.visualization.wcsaxes import WCSAxes

# The next two are not used but are called to register functions with external modules
import sunpy.coordinates
import sunpy.io as io
import sunpy.visualization.colormaps
from sunpy import config
from sunpy.coordinates import HeliographicCarrington, HeliographicStonyhurst, get_earth, sun
from sunpy.coordinates.utils import get_rectangle_coordinates
from sunpy.image.resample import resample as sunpy_image_resample
from sunpy.image.resample import reshape_image_to_4d_superpixel
from sunpy.sun import constants
from sunpy.time import is_time, parse_time
from sunpy.util import MetaDict, expand_list
from sunpy.util.decorators import cached_property_based_on, deprecated
from sunpy.util.exceptions import SunpyDeprecationWarning, SunpyMetadataWarning, SunpyUserWarning
from sunpy.util.functools import seconddispatch
from sunpy.visualization import axis_labels_from_ctype, peek_show, wcsaxes_compat

TIME_FORMAT = config.get("general", "time_format")
PixelPair = namedtuple('PixelPair', 'x y')
SpatialPair = namedtuple('SpatialPair', 'axis1 axis2')

_META_FIX_URL = 'https://docs.sunpy.org/en/stable/code_ref/map.html#fixing-map-metadata'

__all__ = ['GenericMap']


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

    Examples
    --------
    >>> import sunpy.map
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
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
            -128.03072  , -128.03072  ]], dtype=float32)

    >>> aia.spatial_units   # doctest: +REMOTE_DATA
    SpatialPair(axis1=Unit("arcsec"), axis2=Unit("arcsec"))
    >>> aia.peek()   # doctest: +SKIP

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

    This class makes some assumptions about the WCS information contained in
    the meta data. The first and most extensive assumption is that it is
    FITS-like WCS information as defined in the FITS WCS papers.

    Within this scope it also makes some other assumptions.

    * In the case of APIS convention headers where the CROTAi/j arguments are
      provided it assumes that these can be converted to the standard PCi_j
      notation using equations 32 in Thompson (2006).

    * If a CDi_j matrix is provided it is assumed that it can be converted to a
      PCi_j matrix and CDELT keywords as described in
      `Greisen & Calabretta (2002) <https://doi.org/10.1051/0004-6361:20021327>`_

    * The 'standard' FITS keywords that are used by this class are the PCi_j
      matrix and CDELT, along with the other keywords specified in the WCS
      papers. All subclasses of this class must convert their header
      information to this formalism. The CROTA to PCi_j conversion is done in
      this class.

    .. warning::
        This class currently assumes that a header with the CDi_j matrix
        information also includes the CDELT keywords, without these keywords
        this class will not process the WCS.
        Also the rotation_matrix does not work if the CDELT1 and CDELT2
        keywords are exactly equal.
        Also, if a file with more than two dimensions is feed into the class,
        only the first two dimensions (NAXIS1, NAXIS2) will be loaded and the
        rest will be discarded.
    """

    _registry = dict()

    def __init_subclass__(cls, **kwargs):
        """
        An __init_subclass__ hook initializes all of the subclasses of a given class.
        So for each subclass, it will call this block of code on import.
        This replicates some metaclass magic without the need to be aware of metaclasses.
        Here we use this to register each subclass in a dict that has the
        `is_datasource_for` attribute.
        This is then passed into the Map Factory so we can register them.
        """
        super().__init_subclass__(**kwargs)
        if hasattr(cls, 'is_datasource_for'):
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
            warnings.warn("This file contains more than 2 dimensions. "
                          "Data will be truncated to the first two dimensions.", SunpyUserWarning)

        super().__init__(data, meta=MetaDict(header), **kwargs)

        # Correct possibly missing meta keywords
        self._fix_date()
        self._fix_naxis()

        # Setup some attributes
        self._nickname = None
        # These are palceholders for default attributes, which are only set
        # once if their data isn't present in the map metadata.
        self._default_time = None
        self._default_dsun = None
        self._default_carrington_longitude = None
        self._default_heliographic_latitude = None
        self._default_heliographic_longitude = None

        # Validate header
        # TODO: This should be a function of the header, not of the map
        self._validate_meta()
        self._shift = SpatialPair(0 * u.arcsec, 0 * u.arcsec)

        if self.dtype == np.uint8:
            norm = None
        else:
            # Put import here to reduce sunpy.map import time
            from matplotlib import colors
            norm = colors.Normalize()

        # Visualization attributes
        self.plot_settings = {'cmap': 'gray',
                              'norm': norm,
                              'interpolation': 'nearest',
                              'origin': 'lower'
                              }
        if plot_settings:
            self.plot_settings.update(plot_settings)

    def __getitem__(self, key):
        """ This should allow indexing by physical coordinate """
        raise NotImplementedError(
            "The ability to index Map by physical"
            " coordinate is not yet implemented.")

    def _text_summary(self):
        return textwrap.dedent("""\
                   SunPy Map
                   ---------
                   Observatory:\t\t {obs}
                   Instrument:\t\t {inst}
                   Detector:\t\t {det}
                   Measurement:\t\t {meas}
                   Wavelength:\t\t {wave}
                   Observation Date:\t {date}
                   Exposure Time:\t\t {dt:f}
                   Dimension:\t\t {dim}
                   Coordinate System:\t {coord}
                   Scale:\t\t\t {scale}
                   Reference Pixel:\t {refpix}
                   Reference Coord:\t {refcoord}\
                   """).format(obs=self.observatory, inst=self.instrument, det=self.detector,
                               meas=self.measurement, wave=self.wavelength,
                               date=self.date.strftime(TIME_FORMAT),
                               dt=self.exposure_time,
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

    def _repr_html_(self):
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
        # doing set_bad()
        cmap = copy.copy(cm.get_cmap('gray'))
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
        ax.hist(finite_data.ravel(), bins=100, histtype='stepfilled')
        ax.set_facecolor('white')
        ax.semilogy()
        # Explicitly set the power limits for the X axis formatter to avoid text overlaps
        ax.xaxis.get_major_formatter().set_powerlimits((-3, 4))
        ax.set_xlabel('Pixel value')
        ax.set_ylabel('# of pixels')
        ax.set_title('Distribution of pixels with finite values')
        hist_src = _figure_to_base64(fig)

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
                    <td><img src='data:image/png;base64,{hist_src}'/></td>
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
        >>> smap = Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
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
                    <body>{self._repr_html_()}</body>
                </html>"""))
        webbrowser.open_new_tab(url)

    @classmethod
    def _new_instance(cls, data, meta, plot_settings=None, **kwargs):
        """
        Instantiate a new instance of this class using given data.
        This is a shortcut for ``type(self)(data, meta, plot_settings)``.
        """
        return cls(data, meta, plot_settings=plot_settings, **kwargs)

    def _get_lon_lat(self, frame):
        """
        Given a coordinate frame, extract the lon and lat by casting to
        SphericalRepresentation first.
        """
        r = frame.represent_as(UnitSphericalRepresentation)
        return r.lon.to(self.spatial_units[0]), r.lat.to(self.spatial_units[1])

    @property
    def _meta_hash(self):
        return self.meta.item_hash()

    @property
    @cached_property_based_on('_meta_hash')
    def wcs(self):
        """
        The `~astropy.wcs.WCS` property of the map.
        """
        # Construct the WCS based on the FITS header, but don't "do_set" which
        # analyses the FITS header for correctness.
        with warnings.catch_warnings():
            # Ignore warnings we may raise when constructing the fits header about dropped keys.
            warnings.simplefilter("ignore", SunpyUserWarning)
            try:
                w2 = astropy.wcs.WCS(header=self.fits_header, _do_set=False)
            except Exception as e:
                warnings.warn("Unable to treat `.meta` as a FITS header, assuming a simple WCS. "
                              f"The exception raised was:\n{e}")
                w2 = astropy.wcs.WCS(naxis=2)

        # If the FITS header is > 2D pick the first 2 and move on.
        # This will require the FITS header to be valid.
        if w2.naxis > 2:
            # We have to change this or else the WCS doesn't parse properly, even
            # though we don't care about the third dimension. This applies to both
            # EIT and IRIS data, it is here to reduce the chances of silly errors.
            if self.meta.get('cdelt3', None) == 0:
                w2.wcs.cdelt[2] = 1e-10

            w2 = w2.sub([1, 2])

        # Add one to go from zero-based to one-based indexing
        w2.wcs.crpix = u.Quantity(self.reference_pixel) + 1 * u.pix
        # Make these a quantity array to prevent the numpy setting element of
        # array with sequence error.
        w2.wcs.cdelt = u.Quantity(self.scale)
        w2.wcs.crval = u.Quantity([self._reference_longitude, self._reference_latitude])
        w2.wcs.ctype = self.coordinate_system
        w2.wcs.pc = self.rotation_matrix
        w2.wcs.cunit = self.spatial_units
        w2.wcs.dateobs = self.date.isot
        w2.wcs.aux.rsun_ref = self.rsun_meters.to_value(u.m)

        # Astropy WCS does not understand the SOHO default of "solar-x" and
        # "solar-y" ctypes.  This overrides the default assignment and
        # changes it to a ctype that is understood.  See Thompson, 2006, A.&A.,
        # 449, 791.
        if w2.wcs.ctype[0].lower() in ("solar-x", "solar_x"):
            w2.wcs.ctype[0] = 'HPLN-TAN'

        if w2.wcs.ctype[1].lower() in ("solar-y", "solar_y"):
            w2.wcs.ctype[1] = 'HPLT-TAN'

        # Set observer coordinate information
        #
        # Clear all the aux information that was set earlier. This is to avoid
        # issues with maps that store multiple observer coordinate keywords.
        # Note that we have to create a new WCS as it's not possible to modify
        # wcs.wcs.aux in place.
        header = w2.to_header()
        for kw in ['crln_obs', 'dsun_obs', 'hgln_obs', 'hglt_obs']:
            header.pop(kw, None)
        w2 = astropy.wcs.WCS(header)

        # Get observer coord, and transform if needed
        obs_coord = self.observer_coordinate
        if not isinstance(obs_coord.frame, (HeliographicStonyhurst, HeliographicCarrington)):
            obs_coord = obs_coord.transform_to(HeliographicStonyhurst(obstime=self.date))

        sunpy.coordinates.wcs_utils._set_wcs_aux_obs_coord(w2, obs_coord)

        # Validate the WCS here.
        w2.wcs.set()
        return w2

    @property
    def coordinate_frame(self):
        """
        An `astropy.coordinates.BaseFrame` instance created from the coordinate
        information for this Map, or None if the frame cannot be determined.
        """
        try:
            return astropy.wcs.utils.wcs_to_celestial_frame(self.wcs)
        except ValueError as e:
            warnings.warn(f'Could not determine coordinate frame from map metadata.\n{e}',
                          SunpyUserWarning)
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
        axes. See https://wcsaxes.readthedocs.io for more information.
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
    @deprecated(since="2.1", message="Use map.data.size instead", alternative="map.data.size")
    def size(self):
        """
        The number of pixels in the array of the map.
        """
        return u.Quantity(self.data.size, 'pixel')

    @property
    def ndim(self):
        """
        The value of `numpy.ndarray.ndim` of the data array of the map.
        """
        return self.data.ndim

    def std(self, *args, **kwargs):
        """
        Calculate the standard deviation of the data array.
        """
        return self.data.std(*args, **kwargs)

    def mean(self, *args, **kwargs):
        """
        Calculate the mean of the data array.
        """
        return self.data.mean(*args, **kwargs)

    def min(self, *args, **kwargs):
        """
        Calculate the minimum value of the data array.
        """
        return self.data.min(*args, **kwargs)

    def max(self, *args, **kwargs):
        """
        Calculate the maximum value of the data array.
        """
        return self.data.max(*args, **kwargs)

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

        unit = u.Unit(unit_str, format='fits', parse_strict='silent')
        if isinstance(unit, u.UnrecognizedUnit):
            warnings.warn(f'Could not parse unit string "{unit_str}" as a valid FITS unit.\n'
                          f'See {_META_FIX_URL} for how to fix metadata before loading it '
                          'with sunpy.map.Map.\n'
                          'See https://fits.gsfc.nasa.gov/fits_standard.html for'
                          'the FITS unit standards.',
                          SunpyMetadataWarning)
            unit = None
        return unit

# #### Keyword attribute and other attribute definitions #### #

    def _base_name(self):
        """Abstract the shared bit between name and latex_name"""
        return "{nickname} {{measurement}} {date}".format(
            nickname=self.nickname,
            date=parse_time(self.date).strftime(TIME_FORMAT)
        )

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

    @property
    def date(self):
        """
        Image observation time.

        This is taken from the 'DATE-OBS' FITS keyword.
        """
        time = self.meta.get('date-obs', None)
        if time is None:
            if self._default_time is None:
                warnings.warn("Missing metadata for observation time, "
                              "setting observation time to current time. "
                              "Set the 'DATE-OBS' FITS keyword to prevent this warning.",
                              SunpyUserWarning)
                self._default_time = parse_time('now')
            time = self._default_time
        return parse_time(time)

    @property
    def detector(self):
        """
        Detector name.

        This is taken from the 'DETECTOR' FITS keyword.
        """
        return self.meta.get('detector', "")

    @property
    def exposure_time(self):
        """
        Exposure time of the image in seconds.

        This is taken from the 'EXPTIME' FITS keyword.
        """
        return self.meta.get('exptime', 0.0) * u.s

    @property
    def instrument(self):
        """Instrument name."""
        return self.meta.get('instrume', "").replace("_", " ")

    @property
    def measurement(self):
        """
        Measurement name, defaults to the wavelength of image.

        This is taken from the 'WAVELNTH' FITS keyword.
        """
        return u.Quantity(self.meta.get('wavelnth', 0),
                          self.waveunit)

    @property
    def waveunit(self):
        """
        The `~astropy.units.Unit` of the wavelength of this observation.

        This is taken from the 'WAVEUNIT' FITS keyword.
        """
        unit = self.meta.get("waveunit")
        if unit is None:
            return u.one
        return u.Unit(unit)

    @property
    def wavelength(self):
        """
        Wavelength of the observation.

        This is taken from the 'WAVELNTH' FITS keyword.
        """
        return u.Quantity(self.meta.get('wavelnth', 0),
                          self.waveunit)

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
        return self.pixel_to_world(0*u.pix, 0*u.pix)

    @property
    def top_right_coord(self):
        """
        The physical coordinate at the center of the the top right ([-1, -1]) pixel.
        """
        top_right = u.Quantity(self.dimensions) - 1 * u.pix
        return self.pixel_to_world(*top_right)

    @property
    def center(self):
        """
        Return a coordinate object for the center pixel of the array.

        If the array has an even number of pixels in a given dimension,
        the coordinate returned lies on the edge between the two central pixels.
        """
        center = (u.Quantity(self.dimensions) - 1 * u.pix) / 2.
        return self.pixel_to_world(*center)

    @property
    def shifted_value(self):
        """The total shift applied to the reference coordinate by past applications of
        `~sunpy.map.GenericMap.shift`."""
        return self._shift

    @u.quantity_input
    def shift(self, axis1: u.deg, axis2: u.deg):
        """
        Returns a map shifted by a specified amount to, for example, correct
        for a bad map location. These values are applied directly to the
        `~sunpy.map.GenericMap.reference_coordinate`. To check how much shift
        has already been applied see `~sunpy.map.GenericMap.shifted_value`

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
        new_meta['crval1'] = ((self.meta['crval1'] *
                               self.spatial_units[0] + axis1).to(self.spatial_units[0])).value
        new_meta['crval2'] = ((self.meta['crval2'] *
                               self.spatial_units[1] + axis2).to(self.spatial_units[1])).value

        # Create new map with the modification
        new_map = self._new_instance(self.data, new_meta, self.plot_settings)

        new_map._shift = SpatialPair(self.shifted_value[0] + axis1,
                                     self.shifted_value[1] + axis2)

        return new_map

    @property
    def rsun_meters(self):
        """Radius of the sun in meters."""
        return u.Quantity(self.meta.get('rsun_ref', constants.radius), 'meter')

    @property
    def rsun_obs(self):
        """
        Angular radius of the Sun.

        Notes
        -----
        This value is taken the ``'rsun_obs'``, ``'solar_r'``, or ``radius``
        FITS keywords. If none of these keys are present the photospheric limb
        as seen from the observer coordinate is returned.
        """
        rsun_arcseconds = self.meta.get('rsun_obs',
                                        self.meta.get('solar_r',
                                                      self.meta.get('radius',
                                                                    None)))

        if rsun_arcseconds is None:
            warnings.warn("Missing metadata for solar angular radius: assuming photospheric limb "
                          "as seen from observer coordinate.",
                          SunpyUserWarning)
            dsun = self.dsun
            rsun = sun._angular_radius(constants.radius, dsun)
        else:
            rsun = rsun_arcseconds * u.arcsec
        return rsun

    @property
    def coordinate_system(self):
        """
        Coordinate system used for x and y axes (ctype1/2).

        If not present, defaults to (HPLN-TAN, HPLT-TAN), and emits a warning.
        """
        ctype1 = self.meta.get('ctype1', None)
        if ctype1 is None:
            warnings.warn("Missing CTYPE1 from metadata, assuming CTYPE1 is HPLN-TAN",
                          SunpyUserWarning)
            ctype1 = 'HPLN-TAN'

        ctype2 = self.meta.get('ctype2', None)
        if ctype2 is None:
            warnings.warn("Missing CTYPE2 from metadata, assuming CTYPE2 is HPLT-TAN",
                          SunpyUserWarning)
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

    def _remove_existing_observer_location(self):
        """
        Remove all keys that this map might use for observer location.
        """
        all_keys = expand_list([e[0] for e in self._supported_observer_coordinates])
        for key in all_keys:
            self.meta.pop(key)

    @property
    def observer_coordinate(self):
        """
        The Heliographic Stonyhurst Coordinate of the observer.
        """
        missing_meta = {}
        for keys, kwargs in self._supported_observer_coordinates:
            meta_list = [k in self.meta for k in keys]
            if all(meta_list):
                sc = SkyCoord(obstime=self.date, **kwargs)

                # If the observer location is supplied in Carrington coordinates,
                # the coordinate's `observer` attribute should be set to "self"
                if isinstance(sc.frame, HeliographicCarrington):
                    sc.frame._observer = "self"

                return sc.heliographic_stonyhurst
            elif any(meta_list) and not set(keys).isdisjoint(self.meta.keys()):
                if not isinstance(kwargs['frame'], str):
                    kwargs['frame'] = kwargs['frame'].name
                missing_meta[kwargs['frame']] = set(keys).difference(self.meta.keys())

        warning_message = "".join(
            [f"For frame '{frame}' the following metadata is missing: {','.join(keys)}\n" for frame, keys in missing_meta.items()])
        warning_message = "Missing metadata for observer: assuming Earth-based observer.\n" + warning_message
        warnings.warn(warning_message, SunpyMetadataWarning, stacklevel=3)

        return get_earth(self.date)

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
        hgc_frame = HeliographicCarrington(observer=self.observer_coordinate, obstime=self.date)
        return self.observer_coordinate.transform_to(hgc_frame).lat

    @property
    def carrington_longitude(self):
        """Observer Carrington longitude."""
        hgc_frame = HeliographicCarrington(observer=self.observer_coordinate, obstime=self.date)
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
        return PixelPair((self.meta.get('crpix1',
                                        (self.meta.get('naxis1') + 1) / 2.) - 1) * u.pixel,
                         (self.meta.get('crpix2',
                                        (self.meta.get('naxis2') + 1) / 2.) - 1) * u.pixel)

    @property
    def scale(self):
        """
        Image scale along the x and y axes in units/pixel
        (i.e. cdelt1, cdelt2).
        """
        # TODO: Fix this if only CDi_j matrix is provided
        return SpatialPair(self.meta.get('cdelt1', 1.) * self.spatial_units[0] / u.pixel,
                           self.meta.get('cdelt2', 1.) * self.spatial_units[1] / u.pixel)

    @property
    def spatial_units(self):
        """
        Image coordinate units along the x and y axes (i.e. cunit1, cunit2).
        """
        return SpatialPair(u.Unit(self.meta.get('cunit1')),
                           u.Unit(self.meta.get('cunit2')))

    @property
    def rotation_matrix(self):
        r"""
        Matrix describing the rotation required to align solar North with
        the top of the image.

        The order or precendence of FITS keywords which this is taken from is:
        - PC\*_\*
        - CD\*_\*
        - CROTA\*
        """
        if 'PC1_1' in self.meta:
            return np.array([[self.meta['PC1_1'], self.meta['PC1_2']],
                             [self.meta['PC2_1'], self.meta['PC2_2']]])

        elif 'CD1_1' in self.meta:
            cd = np.array([[self.meta['CD1_1'], self.meta['CD1_2']],
                           [self.meta['CD2_1'], self.meta['CD2_2']]])

            cdelt = u.Quantity(self.scale).value

            return cd / cdelt
        else:
            return self._rotation_matrix_from_crota()

    def _rotation_matrix_from_crota(self):
        """
        This method converts the deprecated CROTA FITS kwargs to the new
        PC rotation matrix.

        This method can be overriden if an instruments header does not use this
        conversion.
        """
        lam = self.scale[0] / self.scale[1]
        p = np.deg2rad(self.meta.get('CROTA2', 0))

        return np.array([[np.cos(p), -1 * lam * np.sin(p)],
                         [1/lam * np.sin(p), np.cos(p)]])

    @property
    def fits_header(self):
        """
        A `~astropy.io.fits.Header` representation of the ``meta`` attribute.
        """
        return sunpy.io.fits.header_to_fits(self.meta)

# #### Miscellaneous #### #

    def _fix_date(self):
        # Check commonly used but non-standard FITS keyword for observation
        # time and correct the keyword if we can. Keep updating old one for
        # backwards compatibility.
        if is_time(self.meta.get('date_obs', None)):
            self.meta['date-obs'] = self.meta['date_obs']

    def _fix_naxis(self):
        # If naxis is not specified, get it from the array shape
        if 'naxis1' not in self.meta:
            self.meta['naxis1'] = self.data.shape[1]
        if 'naxis2' not in self.meta:
            self.meta['naxis2'] = self.data.shape[0]
        if 'naxis' not in self.meta:
            self.meta['naxis'] = self.ndim

    def _fix_bitpix(self):
        # Bit-depth
        #
        #   8    Character or unsigned binary integer
        #  16    16-bit twos-complement binary integer
        #  32    32-bit twos-complement binary integer
        # -32    IEEE single precision floating point
        # -64    IEEE double precision floating point
        #
        if 'bitpix' not in self.meta:
            float_fac = -1 if self.dtype.kind == "f" else 1
            self.meta['bitpix'] = float_fac * 8 * self.dtype.itemsize

    def _get_cmap_name(self):
        """Build the default color map name."""
        cmap_string = (self.observatory + self.detector +
                       str(int(self.wavelength.to('angstrom').value)))
        return cmap_string.lower()

    def _validate_meta(self):
        """
        Validates the meta-information associated with a Map.

        This method includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-specific
        validation should be handled in the relevant file in the
        sunpy.map.sources package.

        Allows for default unit assignment for:
            CUNIT1, CUNIT2, WAVEUNIT

        """
        msg = ('Image coordinate units for axis {} not present in metadata.')
        err_message = []
        for i in [1, 2]:
            if self.meta.get(f'cunit{i}') is None:
                err_message.append(msg.format(i, i))

        if err_message:
            err_message.append(
                f'See {_META_FIX_URL} for instructions on how to add missing metadata.')
            raise MapMetaValidationError('\n'.join(err_message))

        for meta_property in ('waveunit', ):
            if (self.meta.get(meta_property) and
                u.Unit(self.meta.get(meta_property),
                       parse_strict='silent').physical_type == 'unknown'):
                warnings.warn(f"Unknown value for {meta_property.upper()}.", SunpyUserWarning)

        if (self.coordinate_system[0].startswith(('SOLX', 'SOLY')) or
                self.coordinate_system[1].startswith(('SOLX', 'SOLY'))):
            warnings.warn("SunPy Map does not support three dimensional data "
                          "and therefore cannot represent heliocentric coordinates. Proceed at your own risk.",
                          SunpyUserWarning)

# #### Data conversion routines #### #

    @staticmethod
    def _check_origin(origin):
        """
        Check origin is valid, and raise a deprecation warning if it's not None.
        """
        if origin is not None:
            warnings.warn('The origin argument is deprecated. If using origin=1, '
                          'manually subtract 1 from your pixels do not pass a value for origin.',
                          SunpyDeprecationWarning)
        if origin not in [None, 0, 1]:
            raise ValueError('origin must be 0 or 1.')

    def world_to_pixel(self, coordinate, origin=None):
        """
        Convert a world (data) coordinate to a pixel coordinate.

        Parameters
        ----------
        coordinate : `~astropy.coordinates.SkyCoord` or `~astropy.coordinates.BaseFrame`
            The coordinate object to convert to pixel coordinates.

        origin : int
            Deprecated.

            Origin of the top-left corner. i.e. count from 0 or 1.
            Normally, origin should be 0 when passing numpy indices, or 1 if
            passing values from FITS header or map attributes.

        Returns
        -------
        x : `~astropy.units.Quantity`
            Pixel coordinate on the CTYPE1 axis.

        y : `~astropy.units.Quantity`
            Pixel coordinate on the CTYPE2 axis.
        """
        self._check_origin(origin)
        x, y = self.wcs.world_to_pixel(coordinate)
        if origin == 1:
            x += 1
            y += 1

        return PixelPair(x * u.pixel, y * u.pixel)

    @u.quantity_input
    def pixel_to_world(self, x: u.pixel, y: u.pixel, origin=None):
        """
        Convert a pixel coordinate to a data (world) coordinate.

        Parameters
        ----------
        x : `~astropy.units.Quantity`
            Pixel coordinate of the CTYPE1 axis. (Normally solar-x).

        y : `~astropy.units.Quantity`
            Pixel coordinate of the CTYPE2 axis. (Normally solar-y).

        origin : int
            Deprecated.

            Origin of the top-left corner. i.e. count from 0 or 1.
            Normally, origin should be 0 when passing numpy indices, or 1 if
            passing values from FITS header or map attributes.

        Returns
        -------
        coord : `astropy.coordinates.SkyCoord`
            A coordinate object representing the output coordinate.
        """
        self._check_origin(origin)
        if origin == 1:
            x = x - 1 * u.pixel
            y = y - 1 * u.pixel

        return self.wcs.pixel_to_world(x, y)

# #### I/O routines #### #

    def save(self, filepath, filetype='auto', **kwargs):
        """Saves the SunPy Map object to a file.

        Currently SunPy can only save files in the FITS format. In the future
        support will be added for saving to other formats.

        Parameters
        ----------
        filepath : str
            Location to save file to.
        filetype : str
            'auto' or any supported file extension.
        hdu_type: None, `~astropy.io.fits.CompImageHDU`
            `None` will return a normal FITS file.
            `~astropy.io.fits.CompImageHDU` will rice compress the FITS file.
        kwargs :
            Any additional keyword arguments are passed to
            `~sunpy.io.write_file`.
        """
        io.write_file(filepath, self.data, self.meta, filetype=filetype,
                      **kwargs)

# #### Image processing routines #### #

    @u.quantity_input
    def resample(self, dimensions: u.pixel, method='linear'):
        """
        Resample to new dimension sizes.

        Uses the same parameters and creates the same co-ordinate lookup points
        as IDL''s congrid routine, which apparently originally came from a
        VAX/VMS routine of the same name.

        Parameters
        ----------
        dimensions : `~astropy.units.Quantity`
            Output pixel dimensions. The first argument corresponds to the 'x'
            axis and the second argument corresponds to the 'y' axis.
        method : str
            Method to use for resampling interpolation.
                * ``'neighbor'`` - Take closest value from original data.
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
        `Rebinning <https://scipy-cookbook.readthedocs.io/items/Rebinning.html>`_
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
        new_meta['cdelt1'] *= scale_factor_x
        new_meta['cdelt2'] *= scale_factor_y
        if 'CD1_1' in new_meta:
            new_meta['CD1_1'] *= scale_factor_x
            new_meta['CD2_1'] *= scale_factor_x
            new_meta['CD1_2'] *= scale_factor_y
            new_meta['CD2_2'] *= scale_factor_y
        new_meta['crpix1'] = (dimensions[0].value + 1) / 2.
        new_meta['crpix2'] = (dimensions[1].value + 1) / 2.
        lon, lat = self._get_lon_lat(self.center.frame)
        new_meta['crval1'] = lon.value
        new_meta['crval2'] = lat.value
        new_meta['naxis1'] = new_data.shape[1]
        new_meta['naxis2'] = new_data.shape[0]

        # Create new map instance
        new_map = self._new_instance(new_data, new_meta, self.plot_settings)
        return new_map

    @u.quantity_input
    def rotate(self, angle: u.deg = None, rmatrix=None, order=4, scale=1.0,
               recenter=False, missing=0.0, use_scipy=False):
        """
        Returns a new rotated and rescaled map.

        Specify either a rotation angle or a rotation matrix, but not both. If
        neither an angle or a rotation matrix are specified, the map will be
        rotated by the rotation angle in the metadata.

        The map will be rotated around the reference coordinate defined in the
        meta data.

        This method also updates the ``rotation_matrix`` attribute and any
        appropriate header data so that they correctly describe the new map.

        Parameters
        ----------
        angle : `~astropy.units.Quantity`
            The angle (degrees) to rotate counterclockwise.
        rmatrix : ndarray
            2x2 linear transformation rotation matrix.
        order : int
            Interpolation order to be used. Must be in the range 0-5.
            When using scikit-image this
            parameter is passed into :func:`skimage.transform.warp` (e.g., 4
            corresponds to bi-quartic interpolation).
            When using scipy it is passed into
            :func:`scipy.ndimage.affine_transform` where it
            controls the order of the spline. Faster performance may be
            obtained at the cost of accuracy by using lower values.
            Default: 4
        scale : float
            A scale factor for the image, default is no scaling
        recenter : bool
            If True, position the axis of rotation at the center of the new map
            Default: False
        missing : float
            The numerical value to fill any missing points after rotation.
            Default: 0.0
        use_scipy : bool
            If True, forces the rotation to use
            :func:`scipy.ndimage.affine_transform`, otherwise it
            uses the :func:`skimage.transform.warp`.
            Default: False, unless scikit-image can't be imported

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new Map instance containing the rotated and rescaled data of the
            original map.

        See Also
        --------
        sunpy.image.transform.affine_transform : The routine this method calls
        for the rotation.

        Notes
        -----
        This function will remove old CROTA keywords from the header.
        This function will also convert a CDi_j matrix to a PCi_j matrix.

        See :func:`sunpy.image.transform.affine_transform` for details on the
        transformations, situations when the underlying data is modified prior
        to rotation, and differences from IDL's rot().
        """
        # Put the import here to reduce sunpy.map import time
        from sunpy.image.transform import affine_transform

        if angle is not None and rmatrix is not None:
            raise ValueError("You cannot specify both an angle and a rotation matrix.")
        elif angle is None and rmatrix is None:
            rmatrix = self.rotation_matrix

        if order not in range(6):
            raise ValueError("Order must be between 0 and 5.")

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

        # Calculate the shape in pixels to contain all of the image data
        extent = np.max(np.abs(np.vstack((self.data.shape @ rmatrix,
                                          self.data.shape @ rmatrix.T))), axis=0)

        # Calculate the needed padding or unpadding
        diff = np.asarray(np.ceil((extent - self.data.shape) / 2), dtype=int).ravel()
        # Pad the image array
        pad_x = int(np.max((diff[1], 0)))
        pad_y = int(np.max((diff[0], 0)))

        new_data = np.pad(self.data,
                          ((pad_y, pad_y), (pad_x, pad_x)),
                          mode='constant',
                          constant_values=(missing, missing))
        new_meta['crpix1'] += pad_x
        new_meta['crpix2'] += pad_y

        # All of the following pixel calculations use a pixel origin of 0

        pixel_array_center = (np.flipud(new_data.shape) - 1) / 2.0

        # Create a temporary map so we can use it for the data to pixel calculation.
        temp_map = self._new_instance(new_data, new_meta, self.plot_settings)

        # Convert the axis of rotation from data coordinates to pixel coordinates
        pixel_rotation_center = u.Quantity(temp_map.world_to_pixel(self.reference_coordinate)).value
        del temp_map

        if recenter:
            pixel_center = pixel_rotation_center
        else:
            pixel_center = pixel_array_center

        # Apply the rotation to the image data
        new_data = affine_transform(new_data.T,
                                    np.asarray(rmatrix),
                                    order=order, scale=scale,
                                    image_center=np.flipud(pixel_center),
                                    recenter=recenter, missing=missing,
                                    use_scipy=use_scipy).T

        if recenter:
            new_reference_pixel = pixel_array_center
        else:
            # Calculate new pixel coordinates for the rotation center
            new_reference_pixel = pixel_center + np.dot(rmatrix,
                                                        pixel_rotation_center - pixel_center)
            new_reference_pixel = np.array(new_reference_pixel).ravel()

        # Define the new reference_pixel
        new_meta['crval1'] = rotation_center[0].value
        new_meta['crval2'] = rotation_center[1].value
        new_meta['crpix1'] = new_reference_pixel[0] + 1  # FITS pixel origin is 1
        new_meta['crpix2'] = new_reference_pixel[1] + 1  # FITS pixel origin is 1

        # Unpad the array if necessary
        unpad_x = -np.min((diff[1], 0))
        if unpad_x > 0:
            new_data = new_data[:, unpad_x:-unpad_x]
            new_meta['crpix1'] -= unpad_x
        unpad_y = -np.min((diff[0], 0))
        if unpad_y > 0:
            new_data = new_data[unpad_y:-unpad_y, :]
            new_meta['crpix2'] -= unpad_y

        # Calculate the new rotation matrix to store in the header by
        # "subtracting" the rotation matrix used in the rotate from the old one
        # That being calculate the dot product of the old header data with the
        # inverse of the rotation matrix.
        pc_C = np.dot(self.rotation_matrix, np.linalg.inv(rmatrix))
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

        Parameters
        ----------
        bottom_left : `astropy.units.Quantity` or `~astropy.coordinates.SkyCoord`
            The bottom-left coordinate of the rectangle. If a `~astropy.coordinates.SkyCoord` it can
            have shape ``(2,)`` and simultaneously define ``top_right``. If specifying
            pixel coordinates it must be given as an `~astropy.units.Quantity`
            object with units of `~astropy.units.pix`.
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
        >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
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
        pixel_corners = u.Quantity(self._parse_submap_input(bottom_left, top_right, width, height)).T

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
    def _parse_submap_coord_input(self, bottom_left, top_right, width, height):
        # Use helper function to get top_right as a SkyCoord
        bottom_left, top_right = get_rectangle_coordinates(bottom_left,
                                                           top_right=top_right,
                                                           width=width,
                                                           height=height)
        frame = bottom_left.frame
        left_lon, bottom_lat = self._get_lon_lat(bottom_left)
        right_lon, top_lat = self._get_lon_lat(top_right)
        corners = SkyCoord([left_lon, left_lon, right_lon, right_lon],
                           [bottom_lat, top_lat, top_lat, bottom_lat],
                           frame=frame)

        return tuple(u.Quantity(self.wcs.world_to_pixel(corners), u.pix).T)

    @u.quantity_input
    def superpixel(self, dimensions: u.pixel, offset: u.pixel = (0, 0)*u.pixel, func=np.sum):
        """Returns a new map consisting of superpixels formed by applying
        'func' to the original map data.

        Parameters
        ----------
        dimensions : tuple
            One superpixel in the new map is equal to (dimension[0],
            dimension[1]) pixels of the original map.
            Note: the first argument corresponds to the 'x' axis and the second
            argument corresponds to the 'y' axis.
        offset : tuple
            Offset from (0,0) in original map pixels used to calculate where
            the data used to make the resulting superpixel map starts.
        func :
            Function applied to the original data.
            The function 'func' must take a numpy array as its first argument,
            and support the axis keyword with the meaning of a numpy axis
            keyword (see the description of `~numpy.sum` for an example.)
            The default value of 'func' is `~numpy.sum`; using this causes
            superpixel to sum over (dimension[0], dimension[1]) pixels of the
            original map.

        Returns
        -------
        out : `~sunpy.map.GenericMap` or subclass
            A new Map which has superpixels of the required size.

        References
        ----------
        | `Summarizing blocks of an array using a moving window <https://mail.scipy.org/pipermail/numpy-discussion/2010-July/051760.html>`_
        """

        # Note: because the underlying ndarray is transposed in sense when
        #   compared to the Map, the ndarray is transposed, resampled, then
        #   transposed back.
        # Note: "center" defaults to True in this function because data
        #   coordinates in a Map are at pixel centers.

        if (offset.value[0] < 0) or (offset.value[1] < 0):
            raise ValueError("Offset is strictly non-negative.")

        # Make a copy of the original data, perform reshaping, and apply the
        # function.
        if self.mask is not None:
            reshaped = reshape_image_to_4d_superpixel(np.ma.array(self.data.copy(), mask=self.mask),
                                                      [dimensions.value[1], dimensions.value[0]],
                                                      [offset.value[1], offset.value[0]])
        else:
            reshaped = reshape_image_to_4d_superpixel(self.data.copy(),
                                                      [dimensions.value[1], dimensions.value[0]],
                                                      [offset.value[1], offset.value[0]])
        new_array = func(func(reshaped, axis=3), axis=1)

        # Update image scale and number of pixels

        # create copy of new meta data
        new_meta = self.meta.copy()

        new_nx = new_array.shape[1]
        new_ny = new_array.shape[0]

        # Update metadata
        new_meta['cdelt1'] = (dimensions[0] * self.scale[0]).value
        new_meta['cdelt2'] = (dimensions[1] * self.scale[1]).value
        if 'CD1_1' in new_meta:
            new_meta['CD1_1'] *= dimensions[0].value
            new_meta['CD2_1'] *= dimensions[0].value
            new_meta['CD1_2'] *= dimensions[1].value
            new_meta['CD2_2'] *= dimensions[1].value
        new_meta['crpix1'] = (new_nx + 1) / 2.
        new_meta['crpix2'] = (new_ny + 1) / 2.
        lon, lat = self._get_lon_lat(self.center.frame)
        new_meta['crval1'] = lon.to(self.spatial_units[0]).value + 0.5 * \
            (offset[0]*self.scale[0]).to(self.spatial_units[0]).value
        new_meta['crval2'] = lat.to(self.spatial_units[1]).value + 0.5 * \
            (offset[1]*self.scale[1]).to(self.spatial_units[1]).value

        # Create new map instance
        if self.mask is not None:
            new_data = np.ma.getdata(new_array)
            new_mask = np.ma.getmask(new_array)
        else:
            new_data = new_array
            new_mask = None

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
    def draw_grid(self, axes=None, grid_spacing: u.deg = 15*u.deg, annotate=True, **kwargs):
        """
        Draws a coordinate overlay on the plot in the Heliographic Stonyhurst
        coordinate system.

        To overlay other coordinate systems see the `WCSAxes Documentation
        <https://docs.astropy.org/en/stable/visualization/wcsaxes/overlaying_coordinate_systems.html>`_

        Parameters
        ----------
        axes: `~matplotlib.axes` or `None`
            Axes to plot limb on, or `None` to use current axes.

        grid_spacing: `~astropy.units.Quantity`
            Spacing for longitude and latitude grid, if length two it specifies
            (lon, lat) spacing.

        annotate : `bool`
            Passing `False` disables the axes labels and the ticks on the top and right axes.

        Returns
        -------
        overlay: `~astropy.visualization.wcsaxes.CoordinatesMap`
            The wcsaxes coordinate overlay instance.

        Notes
        -----
        Keyword arguments are passed onto the `sunpy.visualization.wcsaxes_compat.wcsaxes_heliographic_overlay` function.
        """

        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)
        if not wcsaxes_compat.is_wcsaxes(axes):
            raise TypeError("Overlay grids can only be plotted on WCSAxes plots.")
        return wcsaxes_compat.wcsaxes_heliographic_overlay(axes,
                                                           grid_spacing=grid_spacing,
                                                           annotate=annotate,
                                                           **kwargs)

    def draw_limb(self, axes=None, **kwargs):
        """
        Draws a circle representing the solar limb

        Parameters
        ----------
        axes: `~matplotlib.axes` or None
            Axes to plot limb on or None to use current axes.

        Returns
        -------
        circ: list
            A list containing the `~matplotlib.patches.Circle` object that
            has been added to the axes.

        Notes
        -----
        Keyword arguments are passed onto `matplotlib.patches.Circle`.
        """
        # Put import here to reduce sunpy.map import time
        from matplotlib import patches

        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        transform = wcsaxes_compat.get_world_transform(axes)
        if wcsaxes_compat.is_wcsaxes(axes):
            radius = self.rsun_obs.to(u.deg).value
        else:
            radius = self.rsun_obs.value
        c_kw = {'radius': radius,
                'fill': False,
                'color': 'white',
                'zorder': 100,
                'transform': transform
                }
        c_kw.update(kwargs)

        circ = patches.Circle([0, 0], **c_kw)
        axes.add_artist(circ)

        return [circ]

    @u.quantity_input
    def draw_rectangle(self, bottom_left, *, width: u.deg = None, height: u.deg = None,
                       axes=None, top_right=None, **kwargs):
        """
        Draw a rectangle defined in world coordinates on the plot.

        This draws a rectangle that has corners at ``(bottom_left, top_right)``,
        and has sides parallel to coordinate axes of the map.

        If ``width`` and ``height`` are specified, they are respectively added to the
        longitude and latitude of the ``bottom_left`` coordinate to calculate a
        ``top_right`` coordinate.

        Parameters
        ----------
        bottom_left : `~astropy.coordinates.SkyCoord`
            The bottom-left coordinate of the rectangle. It can
            have shape ``(2,)`` to simultaneously define ``top_right``.
        top_right : `~astropy.coordinates.SkyCoord`
            The top-right coordinate of the rectangle.
        width : `astropy.units.Quantity`, optional
            The width of the rectangle. Required if ``top_right`` is omitted.
        height : `astropy.units.Quantity`
            The height of the rectangle. Required if ``top_right`` is omitted.
        axes : `matplotlib.axes.Axes`
            The axes on which to plot the rectangle, defaults to the current
            axes.

        Returns
        -------
        rect : `list`
            A list containing the `~matplotlib.patches.Rectangle` object, after
            it has been added to ``axes``.

        Notes
        -----
        Extra keyword arguments to this function are passed through to the
        `~matplotlib.patches.Rectangle` instance.
        """
        if isinstance(top_right, u.Quantity) and isinstance(width, u.Quantity):
            # The decorator assigns the first positional arg to top_right and so on.
            height = width
            width = top_right
            top_right = None

        bottom_left, top_right = get_rectangle_coordinates(bottom_left,
                                                           top_right=top_right,
                                                           width=width,
                                                           height=height)

        width = Longitude(top_right.spherical.lon - bottom_left.spherical.lon)
        height = Latitude(top_right.spherical.lat - bottom_left.spherical.lat)

        if not axes:
            axes = plt.gca()
        if wcsaxes_compat.is_wcsaxes(axes):
            axes_unit = u.deg
        else:
            axes_unit = self.spatial_units[0]

        coord = bottom_left.transform_to(self.coordinate_frame)
        bottom_left = u.Quantity((coord.spherical.lon, coord.spherical.lat), unit=axes_unit).value

        width = width.to(axes_unit).value
        height = height.to(axes_unit).value
        kwergs = {'transform': wcsaxes_compat.get_world_transform(axes),
                  'color': 'white',
                  'fill': False}
        kwergs.update(kwargs)
        rect = plt.Rectangle(bottom_left, width, height, **kwergs)
        axes.add_artist(rect)
        return [rect]

    @u.quantity_input
    def draw_contours(self, levels: u.percent, axes=None, **contour_args):
        """
        Draw contours of the data.

        Parameters
        ----------
        levels : `~astropy.units.Quantity`
            A list of numbers indicating the contours to draw. These are given
            as a percentage of the maximum value of the map data.

        axes : `matplotlib.axes.Axes`
            The axes on which to plot the contours. Defaults to the current
            axes.

        Returns
        -------
        cs : `list`
            The `~matplotlib.contour.QuadContourSet` object, after it has been added to
            ``axes``.

        Notes
        -----
        Extra keyword arguments to this function are passed through to the
        `~matplotlib.pyplot.contour` function.
        """
        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        # TODO: allow for use of direct input of contours but requires units of
        # map flux which is not yet implemented

        cs = axes.contour(self.data, 0.01 * levels.to('percent').value * self.data.max(),
                          **contour_args)
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
            figure.colorbar(im, pad=pad)

        if draw_limb:
            self.draw_limb(axes=axes)

        if draw_grid:
            if grid_spacing is None:
                self.draw_grid(axes=axes)
            else:
                self.draw_grid(axes=axes, grid_spacing=grid_spacing)

        return figure

    @u.quantity_input
    def plot(self, annotate=True, axes=None, title=True,
             clip_interval: u.percent = None, **imshow_kwargs):
        """
        Plots the map object using matplotlib, in a method equivalent
        to plt.imshow() using nearest neighbour interpolation.

        Parameters
        ----------
        annotate : `bool`, optional
            If `True`, the data is plotted at its natural scale; with
            title and axis labels.

        axes: `~matplotlib.axes` or None
            If provided the image will be plotted on the given axes. Else the
            current matplotlib axes will be used.

        title : `str`, `bool`, optional
            The plot title. If `True`, uses the default title for this map.

        clip_interval : two-element `~astropy.units.Quantity`, optional
            If provided, the data will be clipped to the percentile interval bounded by the two
            numbers.

        **imshow_kwargs  : `dict`
            Any additional imshow arguments are passed to `~matplotlib.axes.Axes.imshow`.

        Examples
        --------
        >>> # Simple Plot with color bar
        >>> aia.plot()   # doctest: +SKIP
        >>> plt.colorbar()   # doctest: +SKIP

        >>> # Add a limb line and grid
        >>> aia.plot()   # doctest: +SKIP
        >>> aia.draw_limb()   # doctest: +SKIP
        >>> aia.draw_grid()   # doctest: +SKIP

        """
        # Get current axes
        if not axes:
            axes = wcsaxes_compat.gca_wcs(self.wcs)

        if not wcsaxes_compat.is_wcsaxes(axes):
            warnings.warn("WCSAxes not being used as the axes object for this plot."
                          " Plots may have unexpected behaviour. To fix this pass "
                          "'projection=map' when creating the axes",
                          SunpyUserWarning)
            # Check if the image is properly oriented
            if not np.array_equal(self.rotation_matrix, np.identity(2)):
                warnings.warn("The axes of this map are not aligned to the pixel grid. Plot axes may be incorrect.",
                              SunpyUserWarning)

        # Normal plot
        try:
            plot_settings = copy.deepcopy(self.plot_settings)
        except NotImplementedError:
            # MPL dev at the moment does not support deepcopy and this is a workaround.
            plot_settings = self.plot_settings
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

            axes.set_xlabel(axis_labels_from_ctype(self.coordinate_system[0],
                                                   self.spatial_units[0]))
            axes.set_ylabel(axis_labels_from_ctype(self.coordinate_system[1],
                                                   self.spatial_units[1]))

        if not wcsaxes_compat.is_wcsaxes(axes):
            bl = self._get_lon_lat(self.bottom_left_coord)
            tr = self._get_lon_lat(self.top_right_coord)
            x_range = list(u.Quantity([bl[0], tr[0]]).to(self.spatial_units[0]).value)
            y_range = list(u.Quantity([bl[1], tr[1]]).to(self.spatial_units[1]).value)
            imshow_args.update({'extent': x_range + y_range})

        # Take a deep copy here so that a norm in imshow_kwargs doesn't get modified
        # by setting it's vmin and vmax
        try:
            imshow_args.update(copy.deepcopy(imshow_kwargs))
        except NotImplementedError:
            # MPL dev at the moment does not support deepcopy and this is a workaround.
            imshow_args.update(imshow_kwargs)

        if clip_interval is not None:
            if len(clip_interval) == 2:
                clip_percentages = clip_interval.to('%').value
                vmin, vmax = AsymmetricPercentileInterval(*clip_percentages).get_limits(self.data)
            else:
                raise ValueError("Clip percentile interval must be specified as two numbers.")

            imshow_args['vmin'] = vmin
            imshow_args['vmax'] = vmax

        if 'norm' in imshow_args:
            norm = imshow_args['norm']
            if 'vmin' in imshow_args:
                if norm.vmin is not None:
                    raise ValueError('Cannot manually specify vmin, as the norm '
                                     'already has vmin set')
                norm.vmin = imshow_args.pop('vmin')
            if 'vmax' in imshow_args:
                if norm.vmax is not None:
                    raise ValueError('Cannot manually specify vmax, as the norm '
                                     'already has vmax set')
                norm.vmax = imshow_args.pop('vmax')

        if self.mask is None:
            ret = axes.imshow(self.data, **imshow_args)
        else:
            ret = axes.imshow(np.ma.array(np.asarray(self.data), mask=self.mask), **imshow_args)

        if wcsaxes_compat.is_wcsaxes(axes):
            wcsaxes_compat.default_wcs_grid(axes)

        # Set current axes/image if pyplot is being used (makes colorbar work)
        for i in plt.get_fignums():
            if axes in plt.figure(i).axes:
                plt.sca(axes)
                plt.sci(ret)

        return ret

    def contour(self, level, **kwargs):
        """
        Returns coordinates of the contours for a given level value.

        For details of the contouring algorithm see `skimage.measure.find_contours`.

        Parameters
        ----------
        level : float, astropy.units.Quantity
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
        >>> aia = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)  # doctest: +REMOTE_DATA
        >>> contours = aia.contour(50000 * u.ct)  # doctest: +REMOTE_DATA
        >>> print(contours[0])  # doctest: +REMOTE_DATA
            <SkyCoord (Helioprojective: obstime=2011-06-07T06:33:02.770, rsun=696000000.0 m, observer=<HeliographicStonyhurst Coordinate (obstime=2011-06-07T06:33:02.770): (lon, lat, radius) in (deg, deg, m)
        (-0.00406308, 0.04787238, 1.51846026e+11)>): (Tx, Ty) in arcsec
        [(719.59798458, -352.60839064), (717.19243987, -353.75348121),
        ...

        See also
        --------
        `skimage.measure.find_contours`
        """
        from skimage import measure

        if self.unit is not None:
            try:
                level = level.to_value(self.unit)
            # Catch cases where level can't be converted or isn't a Quantity
            except (AttributeError, u.UnitConversionError) as e:
                raise u.UnitsError(f'level must be an astropy quantity convertible to {self.unit}') from e

        contours = measure.find_contours(self.data, level=level, **kwargs)
        contours = [self.wcs.array_index_to_world(c[:, 0], c[:, 1]) for c in contours]
        return contours


class InvalidHeaderInformation(ValueError):
    """Exception to raise when an invalid header tag value is encountered for a
    FITS/JPEG 2000 file."""


def _figure_to_base64(fig):
    # Converts a matplotlib Figure to a base64 UTF-8 string
    buf = BytesIO()
    fig.savefig(buf, format='png', facecolor='none')  # works better than transparent=True
    return b64encode(buf.getvalue()).decode('utf-8')
