"""
TimeSeries is a generic time series class from which all other TimeSeries
classes inherit from.
"""

from __future__ import absolute_import, division, print_function

import os.path
import shutil
import warnings
import inspect
from datetime import datetime
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
import pandas

from sunpy import config
from sunpy.time import is_time, TimeRange, parse_time
from sunpy.util.cond_dispatch import ConditionalDispatch, run_cls
from sunpy.extern.six.moves import urllib
from sunpy.extern import six
from sunpy.sun import sun

import astropy.units as u

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")


# pylint: disable=E1101,E1121,W0404,W0612,W0613
__authors__ = ["Alex Hamilton"]
__email__ = "####"


#def GenericTimeSeries(*args, source=None, concatenate=True, **kwargs):
class GenericTimeSeries:
    """
    A generic time series object.

    Attributes
    ----------
    meta : `str` or `dict`
        The comment string or header associated with the data.
    data : `~pandas.DataFrame`
        An pandas DataFrame prepresenting one or more fields as a function of time.

    Parameters
    ----------------
    filename: string or File
        A file to read.
    source: string
        A string identifier for a registered subclass, matched by that
         subclasses `_is_source_for` method.
    concatenate :  boolean
        Concatenate all files into one Lightcurve object if True, or return
        one Lightcurve for each file if False.

    All other keywords are passed to _is_source_for and then __init__.

    Examples ########
    --------
    >>> from sunpy.lightcurve import LightCurve
    >>> import datetime
    >>> import numpy as np
    >>> base = datetime.datetime.today()
    >>> dates = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    >>> intensity = np.sin(np.arange(0, 12 * np.pi, step=(12 * np.pi) / 24 * 60))
    >>> light_curve = LightCurve.create({"param1": intensity}, index=dates)
    >>> light_curve.peek()   # doctest: +SKIP

    References
    ----------
    * `Pandas Documentation <http://pandas.pydata.org/pandas-docs/dev/dsintro.html>`_

    """

    def __init__(self, data, meta=None, units=None, **kwargs):
        self.data = pandas.DataFrame(data)
        if meta == '' or meta is None:
            self.meta = OrderedDict()
            self.meta.update({'name':None})
        else:
            self.meta = OrderedDict(meta)

        if units == '' or units is None:
            self.units = {}
        else:
            self.units = units

        # Validate input data
        #self._validate_meta()
        self._validate_units()

        # Setup some attributes
        self._nickname = self.detector

        """
        #### From the GenericMap class #####
        super(GenericMap, self).__init__(data, meta=header, **kwargs)

        # Correct possibly missing meta keywords
        self._fix_date()
        self._fix_naxis()

        # Setup some attributes
        self._nickname = self.detector

        # Validate header
        # TODO: This should be a function of the header, not of the map
        self._validate_meta()
        self._validate_units()
        self._shift = Pair(0 * u.arcsec, 0 * u.arcsec)

        if self.dtype == np.uint8:
            norm = None
        else:
            norm = colors.Normalize()
        # Visualization attributes
        self.plot_settings = {'cmap': cm.gray,
                              'norm': norm,
                              'interpolation': 'nearest',
                              'origin': 'lower'
                              }
        #### End from the GenericMap class #####
        """

# #### Keyword attribute and other attribute definitions #### #

    @property
    def name(self):
        """Human-readable description of time-series-type"""
        return "{obs} {detector} {instrument} {measurement} {date:{tmf}}".format(obs=self.observatory,
                                                                detector=self.detector,
                                                                instrument=self.instrument,
                                                                measurement=self.measurement,
                                                                date=parse_time(self.date),
                                                                tmf=TIME_FORMAT)

    @property
    def nickname(self):
        """An abbreviated human-readable description of the time-series-type"""
        return self._nickname

    @nickname.setter
    def nickname(self, n):
        self._nickname = n

    @property
    def date(self):
        """Image observation time"""
        time = parse_time(self.meta.get('date-obs', 'now'))
        if time is None:
            warnings.warn_explicit("Missing metadata for observation time. Using current time.",
                                       Warning, __file__, inspect.currentframe().f_back.f_lineno)
        return parse_time(time)

#    @date.setter
#    def date(self, new_date):
#        self.meta['date-obs'] = new_date
#        #propagate change to malformed FITS keywords
#        if is_time(self.meta.get('date_obs', None)):
#            self.meta['date_obs'] = new_date

    @property
    def detector(self):
        """Detector name"""
        return self.meta.get('detector', "")

    @property
    def dsun(self):
        """The observer distance from the Sun."""
        dsun = self.meta.get('dsun_obs', None)

        if dsun is None:
            warnings.warn_explicit("Missing metadata for Sun-spacecraft separation: assuming Sun-Earth distance",
                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
            dsun = sun.sunearth_distance(self.date).to(u.m)

        return u.Quantity(dsun, 'm')

    @property
    def exposure_time(self):
        """Exposure time of the image in seconds."""
        return self.meta.get('exptime', 0.0) * u.s

    @property
    def instrument(self):
        """Instrument name"""
        return self.meta.get('instrume', "").replace("_", " ")

    @property
    def measurement(self):
        """Measurement name, defaults to the wavelength of image"""
        #return u.Quantity(self.meta.get('wavelnth', 0), self.meta.get('waveunit', ""))
        return 'measurment'

#    @property
#    def wavelength(self):
#        """wavelength of the observation"""
#        return u.Quantity(self.meta.get('wavelnth', 0), self.meta.get('waveunit', ""))

    @property
    def observatory(self):
        """Observatory or Telescope name"""
        return self.meta.get('obsrvtry', self.meta.get('telescop', "")).replace("_", " ")


#    @property
#    def xrange(self):
#        """Return the X range of the image from edge to edge."""
#        #TODO: This should be reading from the WCS object
#        xmin = self.center.x - self.dimensions[0] / 2. * self.scale.x
#        xmax = self.center.x + self.dimensions[0] / 2. * self.scale.x
#        return u.Quantity([xmin, xmax])

#    @property
#    def yrange(self):
#        """Return the Y range of the image from edge to edge."""
#        #TODO: This should be reading from the WCS object
#        ymin = self.center.y - self.dimensions[1] / 2. * self.scale.y
#        ymax = self.center.y + self.dimensions[1] / 2. * self.scale.y
#        return u.Quantity([ymin, ymax])

#    @property
#    def center(self):
#        """The offset between the center of the Sun and the center of the map."""
#        return Pair(wcs.get_center(self.dimensions[0], self.scale.x,
#                                   self.reference_pixel.x,
#                                   self.reference_coordinate.x),
#                    wcs.get_center(self.dimensions[1], self.scale.y,
#                                   self.reference_pixel.y,
#                                   self.reference_coordinate.y))

#    @property
#    def shifted_value(self):
#        """The total shift applied to the reference coordinate by past applications of
#        `~sunpy.map.GenericMap.shift`."""
#        return self._shift

#    @u.quantity_input(x=u.deg, y=u.deg)
#    def shift(self, x, y):
#        """Returns a map shifted by a specified amount to, for example, correct for a bad
#        map location. These values are applied directly to the `~sunpy.map.GenericMap.reference_coordinate`.
#        To check how much shift has already been applied see `~sunpy.map.GenericMap.shifted_value`
#        Parameters
#        ----------
#        x : `~astropy.units.Quantity`
#            The shift to apply to the X coordinate.
#        y : `~astropy.units.Quantity`
#            The shift to apply to the Y coordinate
#        Returns
#        -------
#        out : `~sunpy.map.GenericMap` or subclass
#            A new shifted Map.
#        """
#        new_map = deepcopy(self)
#        new_map._shift = Pair(self.shifted_value.x + x, self.shifted_value.y + y)
#
#        new_meta = self.meta.copy()
#
#        # Update crvals
#        new_meta['crval1'] = ((self.meta['crval1'] * self.spatial_units.x + x).to(self.spatial_units.x)).value
#        new_meta['crval2'] = ((self.meta['crval2'] * self.spatial_units.y + y).to(self.spatial_units.y)).value
#
#        new_map.meta = new_meta
#
#        return new_map

#    @property
#    def rsun_meters(self):
#        """Radius of the sun in meters"""
#        return u.Quantity(self.meta.get('rsun_ref', constants.radius), 'meter')

#    @property
#    def rsun_obs(self):
#        """Radius of the Sun."""
#        rsun_arcseconds = self.meta.get('rsun_obs',
#                                        self.meta.get('solar_r',
#                                                      self.meta.get('radius', None)))
#
#        if rsun_arcseconds is None:
#            warnings.warn_explicit("Missing metadata for solar radius: assuming photospheric limb as seen from Earth",
#                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
#            rsun_arcseconds = sun.solar_semidiameter_angular_size(self.date).to('arcsec').value
#
#        return u.Quantity(rsun_arcseconds, 'arcsec')

#    @property
#    def coordinate_system(self):
#        """Coordinate system used for x and y axes (ctype1/2)"""
#        return Pair(self.meta.get('ctype1', 'HPLN-TAN'),
#                    self.meta.get('ctype2', 'HPLT-TAN'))

#    @property
#    def carrington_longitude(self):
#        """Carrington longitude (crln_obs)"""
#        carrington_longitude = self.meta.get('crln_obs', None)
#
#        if carrington_longitude is None:
#            warnings.warn_explicit("Missing metadata for Carrington longitude: assuming Earth-based observer",
#                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
#            carrington_longitude = (sun.heliographic_solar_center(self.date))[0]
#
#        return u.Quantity(carrington_longitude, 'deg')

#    @property
#    def heliographic_latitude(self):
#        """Heliographic latitude"""
#        heliographic_latitude = self.meta.get('hglt_obs',
#                                              self.meta.get('crlt_obs',
#                                                            self.meta.get('solar_b0', None)))
#
#        if heliographic_latitude is None:
#            warnings.warn_explicit("Missing metadata for heliographic latitude: assuming Earth-based observer",
#                                   Warning, __file__, inspect.currentframe().f_back.f_lineno)
#            heliographic_latitude = (sun.heliographic_solar_center(self.date))[1]
#
#        return u.Quantity(heliographic_latitude, 'deg')

#    @property
#    def heliographic_longitude(self):
#        """Heliographic longitude"""
#        return u.Quantity(self.meta.get('hgln_obs', 0.), 'deg')

#    @property
#    def reference_coordinate(self):
#        """Reference point WCS axes in data units (i.e. crval1, crval2). This value
#        includes a shift if one is set."""
#        return Pair(self.meta.get('crval1', 0.) * self.spatial_units.x,
#                    self.meta.get('crval2', 0.) * self.spatial_units.y)

#    @property
#    def reference_pixel(self):
#        """Reference point axes in pixels (i.e. crpix1, crpix2)"""
#        return Pair(self.meta.get('crpix1', (self.meta.get('naxis1') + 1) / 2.) * u.pixel,
#                    self.meta.get('crpix2', (self.meta.get('naxis2') + 1) / 2.) * u.pixel)

#    @property
#    def scale(self):
#        """Image scale along the x and y axes in units/pixel (i.e. cdelt1, cdelt2)"""
#        #TODO: Fix this if only CDi_j matrix is provided
#        return Pair(self.meta.get('cdelt1', 1.) * self.spatial_units.x / u.pixel,
#                    self.meta.get('cdelt2', 1.) * self.spatial_units.y / u.pixel)

#    @property
#    def spatial_units(self):
#        """Image coordinate units along the x and y axes (i.e. cunit1, cunit2)."""
#        return Pair(u.Unit(self.meta.get('cunit1', 'arcsec')),
#                    u.Unit(self.meta.get('cunit2', 'arcsec')))

#    @property
#    def rotation_matrix(self):
#        """Matrix describing the rotation required to align solar North with
#        the top of the image."""
#        if 'PC1_1' in self.meta:
#            return np.matrix([[self.meta['PC1_1'], self.meta['PC1_2']],
#                              [self.meta['PC2_1'], self.meta['PC2_2']]])
#
#        elif 'CD1_1' in self.meta:
#            cd = np.matrix([[self.meta['CD1_1'], self.meta['CD1_2']],
#                            [self.meta['CD2_1'], self.meta['CD2_2']]])
#
#            cdelt = u.Quantity(self.scale).value
#
#            return cd / cdelt
#        else:
#            return self._rotation_matrix_from_crota()

#    def _rotation_matrix_from_crota(self):
#        """
#        This method converts the deprecated CROTA FITS kwargs to the new
#        PC rotation matrix.
#        This method can be overriden if an instruments header does not use this
#        conversion.
#        """
#        lam = self.scale.y / self.scale.x
#        p = np.deg2rad(self.meta.get('CROTA2', 0))
#
#        return np.matrix([[np.cos(p), -1 * lam * np.sin(p)],
#                          [1/lam * np.sin(p), np.cos(p)]])

# #### From Pandas #### #
    def sort_index(self, **kwargs):
        """Returns a truncated version of the lightcurve object.

        Parameters
        ----------
        a : `sunpy.time.TimeRange`
            A time range to truncate to.

        Returns
        -------
        newlc : `~sunpy.lightcurve.LightCurve`
            A new lightcurve with only the selected times.
        """
        return GenericTimeSeries(self.data.sort_index(**kwargs), self.meta.copy(), self.units)

# #### From Generic Spec #### #


# #### From OLD LightCurve #### #
    @property
    def time_range(self):
        """Returns the start and end times of the LightCurve as a `~sunpy.time.TimeRange`
        object"""
        return TimeRange(self.data.index[0], self.data.index[-1])

    def truncate(self, a, b=None, int=None):
        """Returns a truncated version of the lightcurve object.

        Parameters
        ----------
        a : `sunpy.time.TimeRange` or string or ...
            Either a time range to truncate to, or a start time in some format
            recognised by pandas.

        b : string or ...
            If specified, the end time of the time range in some format
            recognised by pandas.

        int : integer
            If specified, the interger indicating the slicing intervals.

        Returns
        -------
        newlc : `~sunpy.lightcurve.LightCurve`
            A new lightcurve with only the selected times.
        """
        if isinstance(a, TimeRange):
            # If we have a TimeRange, extract the values
            start = a.start
            end   = a.end
        else:
            # Otherwise we already have the values
            start = a
            end   = b

        # If an interval integer was given then use in truncation.
        if int:
            truncated = self.data.sort_index()[start:end:int]
        else:
            truncated = self.data.sort_index()[start:end]
        ####return self.__class__.create(truncated, self.meta.copy())
        return GenericTimeSeries(truncated, self.meta.copy(), self.units)

    def extract(self, column_name):
        """Returns a new lightcurve with the chosen column.

        Parameters
        ----------
        column_name : `str`
            A valid column name

        Returns
        -------
        newlc : `~sunpy.lightcurve.LightCurve`
            A new lightcurve with only the selected column.
        """
        """
        # TODO allow the extract function to pick more than one column
        if isinstance(self, pandas.Series):
            return self
        else:
            return LightCurve(self.data[column_name], self.meta.copy())
        """
        # Extract column and remove empty rows
        data = self.data[column_name].dropna()

        # Return a new GenericTimeSeries
        return GenericTimeSeries(data, self.meta.copy(), { column_name:self.units[column_name] })

    def concatenate(self, otherts):
        """Concatenate another light curve. This function will check and remove
        any duplicate times. It will keep the column values from the original
        lightcurve to which the new lightcurve is being added.

        Parameters
        ----------
        otherts : `~sunpy.lightcurve.LightCurve`
            Another time series of the same type.

        Returns
        -------
        newts : `~sunpy.lightcurve.LightCurve`
            A new time series.
        """
        if not isinstance(otherts, self.__class__):
            raise TypeError("Lightcurve classes must match.")

        # Metadata will be the original time series metadata but with an additional
        # entry containing all of the additional time series.
        """# From original implmentation
        meta = OrderedDict()
        meta.update({str(self.data.index[0]):self.meta.copy()})
        meta.update({str(otherlightcurve.data.index[0]):otherlightcurve.meta.copy()})
        """
        meta = self.meta.copy()
        meta['2nd_ts_meta'] = otherts.meta.copy()

        data = pd.concat([self.data.copy(), otherts.data])

        units = OrderedDict()
        units.update(self.units)
        units.update(otherts.units)
        """
        data['index'] = data.index
        # default behavior of drop_duplicates is keep the first column.
        data = data.drop_duplicates(subset='index')
        data.set_index = data['index']
        data.drop('index', axis=1, inplace=True)
        return self.__class__.create(data, meta)
        """
        return GenericTimeSeries(data, meta, units)

# #### Miscellaneous #### #
    def _validate_meta(self):
        """
        Validates the meta-information associated with a TimeSeries.

        This method includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-specific
        validation should be handled in the relevant file in the
        sunpy.timeseries.sources package.

        Allows for default unit assignment for:
            COL_UNITS

        """

        warnings.simplefilter('always', Warning)

        for meta_property in ('cunit1', 'cunit2', 'waveunit'):
            if (self.meta.get(meta_property) and
                u.Unit(self.meta.get(meta_property),
                       parse_strict='silent').physical_type == 'unknown'):

                warnings.warn("Unknown value for "+meta_property.upper(), Warning)

    def _validate_units(self, **kwargs):
        """
        Validates the meta-information associated with a TimeSeries.

        This method includes very basic validation checks which apply to
        all of the kinds of files that SunPy can read. Datasource-specific
        validation should be handled in the relevant file in the
        sunpy.timeseries.sources package.

        Allows for default unit assignment for:
            COL_UNITS

        """

        warnings.simplefilter('always', Warning)

        # For all columns not present in the units dictionary.
        for column in set(self.data.columns.tolist()) - set(self.units.keys()):
            self.units[column] = u.Quantity(1.0)
            warnings.warn("Unknown units for "+column, Warning)

if __name__ == "__main__":
    # Build a traditional lightcurve
    intensity = np.sin(np.arange(0, 12 * np.pi, ((12 * np.pi) / (24*60))))
    import sunpy
    import datetime
    #from sunpy import lightcurve
    dates = [datetime.datetime(2016, 6, 7, 12, 0) - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
    light_curve = sunpy.lightcurve.LightCurve.create({"param1": intensity}, index=dates)

    # Use this to
    units_list = { "param1": u.Quantity(1 * u.meter) }
    gts = GenericTimeSeries(light_curve.data, light_curve.meta, units_list)

    ###########################################################################
    # The basic functions necessary for the API
    ###########################################################################
    # Sorting: sort the data so time is monotonically increasing
    gts.data = gts.data.sort_index()

    # Truncation: getting shorter duration lightcurves from an existing timeseries. (exclusive of end)
    trunc_str_dates = gts.data.sort_index()['20150419':'20160607 01'] # Note: need to sort first

    import datetime
    st = datetime.datetime(2016, 6, 7, 0, 0)
    en = datetime.datetime(2016, 6, 8, 12, 0)
    trunc_datetimes = gts.data.sort_index()[st:en]

    trunc_slice = gts.data[10:500]

    # Subsampling: taking specific elements of an existing lightcurve (for example every n'th element), and creating a new lightcurve
    sub1 = gts.data[::3]
    sub2 = gts.data[1::3]
    sub3 = gts.data[2::3]

    # Merging/Concatenating: merge 2 or more timeseries of the same data
    import pandas as pd
    merged = pd.concat([sub1, sub2, sub3])
    # Will work on different data, as long as the column (title) is different
    sub2.columns = ['param2']
    merged = pd.concat([sub1, sub2, sub3]).sort_index()

    # Resampling: taking an existing lightcurve and resample it at different times to get a new lightcurve
    # summing every 'n' elements of the original time-series
    resampled = gts.data.resample('3T').sum() # '3T = 3 mins'
    resampled = gts.data.resample('60S').mean() # '60S = 60 seconds'
    resampled = gts.data.resample('H').std() # 'H = one hour'
    resampled = gts.data.resample('D').sum() # 'D = one day'
    # Note: use of the methods: .sum(), .mean(), .std()
    # Note: not sure how to do every nth element.

    # Upsampling: increasing the cadence with interpolation.
    ###################upsampled = gts.data[0:5].resample('30S').asfreq() #select first 5 rows
    upsampled = gts.data[0:10].resample('30S').pad()
    upsampled = gts.data[0:10].resample('30S').bfill()
    upsampled = gts.data[0:10].resample('30S').ffill()

    # Unit aware: Being able to assign a physical unit to a column of the LightCurve data
    # Note: only a bolt on ATM, use a self.units ordered dictionary to store quantity related to each column.
    column = 'param1'
    quantity = u.Quantity(gts.data['param1'].values * gts.units[column])


    ###########################################################################
    # Now, the same fuinctions, but acting on the GenericTimeSeries
    ###########################################################################
    # Sorting: sort the data so time is monotonically increasing
    gts = gts.sort_index()
    # Note: not sure the difference here...

    # Truncation: getting shorter duration lightcurves from an existing timeseries. (exclusive of end)
    trunc_str_dates = gts.truncate('20150419','20160607 01')

    tr = TimeRange('2015/04/19', '2016/06/07')
    trunc_time_range = gts.truncate(tr)

    import datetime
    st = datetime.datetime(2016, 6, 7, 0, 0)
    en = datetime.datetime(2016, 6, 8, 12, 0)
    trunc_datetimes = gts.truncate(st,en)

    trunc_slice = gts.truncate(10,500)

    # Subsampling: taking specific elements of an existing lightcurve (for example every n'th element), and creating a new lightcurve
    sub1 = gts.truncate(0,5000,3)
    sub2 = gts.truncate(1,5000,3)
    sub3 = gts.truncate(2,5000,3)

    # Merging/Concatenating: merge 2 or more timeseries of the same data
    import pandas as pd
    merged = sub1.concatenate(sub2)
    # Will work on different data, as long as the column (title) is different
    sub2.data.columns = ['param2']
    merged = sub1.concatenate(sub2).sort_index()

    # Resampling: taking an existing lightcurve and resample it at different times to get a new lightcurve
    # summing every 'n' elements of the original time-series
    resampled = gts.data.resample('3T').sum() # '3T = 3 mins'
    resampled = gts.data.resample('60S').mean() # '60S = 60 seconds'
    resampled = gts.data.resample('H').std() # 'H = one hour'
    resampled = gts.data.resample('D').sum() # 'D = one day'
    # Note: use of the methods: .sum(), .mean(), .std()
    # Note: not sure how to do every nth element.

    # Upsampling: increasing the cadence with interpolation.
    ###################upsampled = gts.data[0:5].resample('30S').asfreq() #select first 5 rows
    upsampled = gts.data[0:10].resample('30S').pad()
    upsampled = gts.data[0:10].resample('30S').bfill()
    upsampled = gts.data[0:10].resample('30S').ffill()

    # Unit aware: Being able to assign a physical unit to a column of the LightCurve data
    # Note: only a bolt on ATM, use a self.units ordered dictionary to store quantity related to each column.
    column = 'param1'
    quantity = u.Quantity(gts.data['param1'].values * gts.units[column])
