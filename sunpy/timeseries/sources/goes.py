"""
This module provies GOES XRS `~sunpy.timeseries.TimeSeries` source.
"""
import datetime
from pathlib import Path
from collections import OrderedDict

import h5netcdf
import matplotlib.dates
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import pyplot as plt
from pandas import DataFrame

import astropy.units as u
from astropy.time import Time, TimeDelta

import sunpy.io
from sunpy.io.file_tools import UnrecognizedFileTypeError
from sunpy.time import is_time_in_given_format, parse_time
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util.metadata import MetaDict
from sunpy.visualization import peek_show

__all__ = ['XRSTimeSeries']


class XRSTimeSeries(GenericTimeSeries):
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

        fig, ax = plt.subplots()

        dates = matplotlib.dates.date2num(parse_time(self.to_dataframe().index).datetime)

        ax.plot_date(dates, self.to_dataframe()['xrsa'], '-',
                     label=r'0.5--4.0 $\AA$', color='blue', lw=2, **kwargs)
        ax.plot_date(dates, self.to_dataframe()['xrsb'], '-',
                     label=r'1.0--8.0 $\AA$', color='red', lw=2, **kwargs)

        ax.set_yscale("log")
        ax.set_ylim(1e-9, 1e-2)
        ax.set_title(title)
        ax.set_ylabel('Watts m$^{-2}$')
        ax.set_xlabel(datetime.datetime.isoformat(self.to_dataframe().index[0])[0:10])

        ax2 = ax.twinx()
        ax2.set_yscale("log")
        ax2.set_ylim(1e-9, 1e-2)
        labels = ['A', 'B', 'C', 'M', 'X']
        centers = np.logspace(-7.5, -3.5, len(labels))
        ax2.yaxis.set_minor_locator(mticker.FixedLocator(centers))
        ax2.set_yticklabels(labels, minor=True)
        ax2.set_yticklabels([])

        ax.yaxis.grid(True, 'major')
        ax.xaxis.grid(False, 'major')
        ax.legend()

        # TODO: display better tick labels for date range (e.g. 06/01 - 06/05)
        formatter = matplotlib.dates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(formatter)

        ax.fmt_xdata = matplotlib.dates.DateFormatter('%H:%M')
        fig.autofmt_xdate()

        return fig

    @classmethod
    def _parse_file(cls, filepath):
        """
        Parses a GOES/XRS FITS file.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        if sunpy.io.detect_filetype(filepath) == "hdf5":
            return cls._parse_netcdf(filepath)
        try:
            hdus = sunpy.io.read_file(filepath)
        except UnrecognizedFileTypeError:
            raise ValueError(f"{Path(filepath).name} is not supported. Only fits and netCDF (nc) can be read.")
        else:
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

        header = MetaDict(OrderedDict(hdulist[0].header))
        if len(hdulist) == 4:
            if is_time_in_given_format(hdulist[0].header['DATE-OBS'], '%d/%m/%Y'):
                start_time = Time.strptime(hdulist[0].header['DATE-OBS'], '%d/%m/%Y')
            elif is_time_in_given_format(hdulist[0].header['DATE-OBS'], '%d/%m/%y'):
                start_time = Time.strptime(hdulist[0].header['DATE-OBS'], '%d/%m/%y')
            else:
                raise ValueError("Date not recognized")
            xrsb = hdulist[2].data['FLUX'][0][:, 0]
            xrsa = hdulist[2].data['FLUX'][0][:, 1]
            seconds_from_start = hdulist[2].data['TIME'][0]
        elif 1 <= len(hdulist) <= 3:
            start_time = parse_time(header['TIMEZERO'], format='utime')
            seconds_from_start = hdulist[0].data[0]
            xrsb = hdulist[0].data[1]
            xrsa = hdulist[0].data[2]
        else:
            raise ValueError("Don't know how to parse this file")

        times = start_time + TimeDelta(seconds_from_start*u.second)
        times.precision = 9

        # remove bad values as defined in header comments
        xrsb[xrsb == -99999] = np.nan
        xrsa[xrsa == -99999] = np.nan

        # fix byte ordering
        newxrsa = xrsa.byteswap().newbyteorder()
        newxrsb = xrsb.byteswap().newbyteorder()

        data = DataFrame({'xrsa': newxrsa, 'xrsb': newxrsb},
                         index=times.isot.astype('datetime64'))
        data.sort_index(inplace=True)

        # Add the units
        units = OrderedDict([('xrsa', u.W/u.m**2),
                             ('xrsb', u.W/u.m**2)])
        return data, header, units

    @staticmethod
    def _parse_netcdf(filepath):
        """
        Parses the netCDF GOES files to return the data, header and associated units.

        Parameters
        ----------
        filepath : `~str`
            The path of the file to parse
        """

        with h5netcdf.File(filepath, mode="r") as d:

            header =  MetaDict(OrderedDict(d.attrs))
            if "a_flux" in d.variables:
                xrsa = np.array(d["a_flux"])
                xrsb = np.array(d["b_flux"])
                start_time_str = d["time"].attrs["units"].astype(str).lstrip("seconds since").rstrip("UTC")
                times = parse_time(start_time_str) + TimeDelta(d["time"], format="sec")
            elif "xrsa_flux" in d.variables:
                xrsa = np.array(d["xrsa_flux"])
                xrsb = np.array(d["xrsb_flux"])
                start_time_str = d["time"].attrs["units"].astype(str).lstrip("seconds since")
                times = parse_time(start_time_str) + TimeDelta(d["time"], format="sec")

            else:
                raise ValueError(f"The file {filepath} doesn't seem to be a GOES netcdf file.")

        data = DataFrame({"xrsa": xrsa, "xrsb": xrsb}, index=times.datetime)
        data = data.replace(-9999, np.nan)
        units = OrderedDict([("xrsa", u.W/u.m**2),
                             ("xrsb", u.W/u.m**2)])

        return data, header, units

    @classmethod
    def is_datasource_for(cls, **kwargs):
        """
        Determines if header corresponds to a GOES lightcurve
        `~sunpy.timeseries.TimeSeries`.
        """
        if "source" in kwargs.keys():
            return kwargs["source"].lower().startswith(cls._source)
        if "meta" in kwargs.keys():
            return kwargs["meta"].get("TELESCOP", "").startswith("GOES")

        if "filepath" in kwargs.keys():
            try:
                if sunpy.io.detect_filetype(kwargs["filepath"]) == "hdf5":
                    with h5netcdf.File(kwargs["filepath"], mode="r") as f:
                        return "XRS" in f.attrs["summary"].astype("str")
            except Exception:
                return False
