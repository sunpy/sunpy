"""
This module provies GOES XRS `~sunpy.timeseries.TimeSeries` source.
"""
import datetime
from collections import OrderedDict

import matplotlib.dates
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import pyplot as plt
from pandas import DataFrame

import astropy.units as u
from astropy.time import Time, TimeDelta

import sunpy.io
from sunpy.time import TimeRange, is_time_in_given_format, parse_time
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util.metadata import MetaDict
from sunpy.visualization import peek_show

__all__ = ['XRSTimeSeries']


class XRSTimeSeries2(GenericTimeSeries):
    """
    GOES XRS Time Series.

    Each GOES satellite there are two X-ray Sensors (XRS) which provide solar X ray fluxes
    for the wavelength bands of 0.5 to 4 Å (short channel) sand 1 to 8 Å (long channel).
    Most recent data is usually available one or two days late.

    Data is available starting on 1981/01/01.

    Examples
    --------
    >>> import sunpy.timeseries
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> goes = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES)  # doctest: +REMOTE_DATA
    >>> goes.peek()   # doctest: +SKIP

    References
    ----------
    * `GOES Mission Homepage <https://www.goes.noaa.gov>`_
    * `GOES XRS Homepage <https://www.swpc.noaa.gov/products/goes-x-ray-flux>`_
    * `GOES XRS Guide <https://ngdc.noaa.gov/stp/satellite/goes/doc/GOES_XRS_readme.pdf>`_
    * `NASCOM Data Archive <https://umbra.nascom.nasa.gov/goes/fits/>`_

    Notes
    -----
    * https://umbra.nascom.nasa.gov/goes/fits/goes_fits_files_notes.txt
    """
    # Class attribute used to specify the source class of the TimeSeries.
    _source = 'xrs'

    @peek_show
    def peek(self, title="GOES Xray Flux", **kwargs):
        """
        Plots GOES XRS light curve is the usual manner. An example is shown
        below:

        .. plot::

            import sunpy.timeseries
            import sunpy.data.sample
            ts_goes = sunpy.timeseries.TimeSeries(sunpy.data.sample.GOES_XRS_TIMESERIES, source='XRS')
            ts_goes.peek()

        Parameters
        ----------
        title : `str`. optional
            The title of the plot. Defaults to "GOES Xray Flux".
        **kwargs : `dict`
            Additional plot keyword arguments that are handed to `axes.plot` functions
        """
        # Check we have a timeseries valid for plotting
        self._validate_data_for_plotting()

        figure = plt.figure()
        axes = plt.gca()

        dates = matplotlib.dates.date2num(parse_time(self.to_dataframe().index).datetime)

        axes.plot_date(dates, self.to_dataframe()['xrsa'], '-',
                       label=r'0.5--4.0 $\AA$', color='blue', lw=2, **kwargs)
        axes.plot_date(dates, self.to_dataframe()['xrsb'], '-',
                       label=r'1.0--8.0 $\AA$', color='red', lw=2, **kwargs)

        axes.set_yscale("log")
        axes.set_ylim(1e-9, 1e-2)
        axes.set_title(title)
        axes.set_ylabel('Watts m$^{-2}$')
        axes.set_xlabel(datetime.datetime.isoformat(self.to_dataframe().index[0])[0:10])

        ax2 = axes.twinx()
        ax2.set_yscale("log")
        ax2.set_ylim(1e-9, 1e-2)
        labels = ['A', 'B', 'C', 'M', 'X']
        centers = np.logspace(-7.5, -3.5, len(labels))
        ax2.yaxis.set_minor_locator(mticker.FixedLocator(centers))
        ax2.set_yticklabels(labels, minor=True)
        ax2.set_yticklabels([])

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(False, 'major')
        axes.legend()

        # TODO: display better tick labels for date range (e.g. 06/01 - 06/05)
        formatter = matplotlib.dates.DateFormatter('%H:%M')
        axes.xaxis.set_major_formatter(formatter)

        axes.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M')
        figure.autofmt_xdate()

        return figure


    @classmethod
    def _parse_file(cls, filepath):
        import xarray
        """
        Parses a GOES/XRS FITS file.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        hdus = xarray.open_dataset(filepath)
        return cls._parse_hdus(hdus)

    @classmethod
    def _parse_hdus(cls, hdulist):
        """
        Parses a GOES/XRS FITS `~astropy.io.fits.HDUList` from a FITS file.

        Parameters
        ----------
        hdulist : `astropy.io.fits.HDUList`
            A HDU list.
        """

        if 'xrsa_flux' in hdulist.keys():
            xrsb = hdulist['xrsa_flux']
            xrsa = hdulist['xrsb_flux']    

        else:  
            xrsb = hdulist['a_flux']
            xrsa = hdulist['b_flux']

        times = parse_time(hdulist['time']).datetime

        data = DataFrame({'xrsa': xrsa, 'xrsb': xrsb},
                         index=times)
        data.sort_index(inplace=True)

        # Add the units
        units = OrderedDict([('xrsa', u.W/u.m**2),
                             ('xrsb', u.W/u.m**2)])

        header = MetaDict(OrderedDict(hdulist.variables))  
        return data, header, units

    @classmethod
    def is_datasource_for(cls, **kwargs):
        """
        Determines if header corresponds to a GOES lightcurve
        `~sunpy.timeseries.TimeSeries`.
        """
        if 'source' in kwargs.keys():
            if kwargs.get('source', ''):
                return kwargs.get('source', '').lower().startswith(cls._source)
        if 'meta' in kwargs.keys():
            return kwargs['meta'].get('TELESCOP', '').startswith('GOES')
