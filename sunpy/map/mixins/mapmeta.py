import re
from collections import namedtuple

import numpy as np

import astropy.units as u
import astropy.wcs
from astropy.coordinates import SkyCoord

import sunpy
from sunpy import config, log
from sunpy.coordinates import HeliographicCarrington, get_earth, sun
from sunpy.sun import constants
from sunpy.time import is_time, parse_time
from sunpy.util import expand_list
from sunpy.util.decorators import cached_property_based_on
from sunpy.util.exceptions import warn_deprecated, warn_metadata, warn_user

__all__ = ['MapMetaValidationError', 'MapMetaMixin', 'PixelPair', 'SpatialPair']


TIME_FORMAT = config.get("general", "time_format")
_META_FIX_URL = 'https://docs.sunpy.org/en/stable/code_ref/map.html#fixing-map-metadata'

PixelPair = namedtuple('PixelPair', 'x y')
SpatialPair = namedtuple('SpatialPair', 'axis1 axis2')


class MapMetaValidationError(AttributeError):
    pass


class MapMetaMixin:
    """
    This class contains the properties for processing metadata
    """

    @property
    def _meta_hash(self):
        return self.meta.item_hash()

    @property
    def coordinate_frame(self):
        """
        An `astropy.coordinates.BaseCoordinateFrame` instance created from the coordinate
        information for this Map, or None if the frame cannot be determined.
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
            wunit = sunpy.io._fits.extract_waveunit(self.meta)
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
        top_right = np.array([self.shape[1], self.shape[0]]) - 1
        return self.wcs.pixel_to_world(*top_right)

    @property
    def center(self):
        """
        Return a coordinate object for the center pixel of the array.

        If the array has an even number of pixels in a given dimension,
        the coordinate returned lies on the edge between the two central pixels.
        """
        center = (np.array([self.shape[1], self.shape[0]]) - 1) / 2.
        return self.wcs.pixel_to_world(*center)

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
            log.debug('Deriving rotation matrix from PC matrix')
            return np.array(
                [
                    [self.meta.get('PC1_1', 1), self.meta.get('PC1_2', 0)],
                    [self.meta.get('PC2_1', 0), self.meta.get('PC2_2', 1)]
                ]
            )
        elif any(key in self.meta for key in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']):
            log.debug('Deriving rotation matrix from CD matrix')
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
            log.debug('Deriving rotation matrix from CROTA')
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
    def _cd_matrix(self):
        """
        Return a CD matrix if the necessary keys exist
        """
        if 'cd1_1' not in self.meta:
            return None
        cd_matrix = np.eye(2)
        for i in range(2):
            for j in range(2):
                cd_matrix[i,j] = self.meta.get(f'CD{i+1}_{j+1}')
        return cd_matrix

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

    @property
    def fits_header(self):
        """
        A `~astropy.io.fits.Header` representation of the ``meta`` attribute.
        """
        return sunpy.io._fits.header_to_fits(self.meta)

    def _validate_meta(self):
        """
        Validates some meta-information associated with a Map.

        This method includes very basic validation checks which apply to
        all of the kinds of files that sunpy can read. Datasource-specific
        validation should be handled in the relevant file in the
        sunpy.map.sources package.
        """
        # Only run this once
        if self._metadata_validated:
            return

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

        self._metadata_validated = True
