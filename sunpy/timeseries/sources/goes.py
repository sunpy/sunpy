"""
This module provides GOES XRS `~sunpy.timeseries.TimeSeries` source.
"""
from pathlib import Path
from collections import OrderedDict

import h5netcdf
import numpy as np
from pandas import DataFrame

import astropy.units as u
from astropy.time import Time, TimeDelta

import sunpy.io
from sunpy import log
from sunpy.extern import parse
from sunpy.io._file_tools import UnrecognizedFileTypeError
from sunpy.time import is_time_in_given_format, parse_time
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util.exceptions import warn_user
from sunpy.util.metadata import MetaDict

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
    * `GOES Mission Homepage <https://www.goes.noaa.gov>`__
    * `GOES XRS Homepage <https://www.swpc.noaa.gov/products/goes-x-ray-flux>`__
    * `GOES XRS Guide <https://ngdc.noaa.gov/stp/satellite/goes/doc/GOES_XRS_readme.pdf>`__
    * `NASCOM Data Archive <https://umbra.nascom.nasa.gov/goes/fits/>`__

    Notes
    -----
    * https://umbra.nascom.nasa.gov/goes/fits/goes_fits_files_notes.txt
    """
    # Class attributes used to specify the source class of the TimeSeries
    # and a URL to the mission website.
    _source = 'xrs'
    _url = "https://www.swpc.noaa.gov/products/goes-x-ray-flux"

    _peek_title = "GOES X-ray flux"

    _netcdf_read_kw = {}
    _netcdf_read_kw['decode_vlen_strings'] = True

    def plot(self, axes=None, columns=None, **kwargs):
        """
        Plots the GOES XRS light curve.

        Parameters
        ----------
        axes : `matplotlib.axes.Axes`, optional
            The axes on which to plot the TimeSeries. Defaults to current axes.
        columns : list[str], optional
            If provided, only plot the specified columns.
        **kwargs : `dict`
            Additional plot keyword arguments that are handed to `~matplotlib.axes.Axes.plot`
            functions.

        Returns
        -------
        `~matplotlib.axes.Axes`
            The plot axes.
        """
        if columns is None:
            columns = ["xrsa", "xrsb"]
        axes, columns = self._setup_axes_columns(axes, columns)
        plot_settings = {"xrsa": ["blue", r"0.5$-$4.0 $\mathrm{\AA}$"], "xrsb": ["red", r"1.0$-$8.0 $\mathrm{\AA}$"]}
        data = self.to_dataframe()
        for channel in columns:
            axes.plot(
                data.index, data[channel], "-", label=plot_settings[channel][1], color=plot_settings[channel][0], lw=2, **kwargs
            )
        axes.set_yscale("log")
        axes.set_ylim(1e-9, 1e-2)
        axes.set_ylabel("Watts m$^{-2}$")
        self._setup_x_axis(axes)

        labels = ['A', 'B', 'C', 'M', 'X']
        centers = np.logspace(-7.5, -3.5, len(labels))
        for value, label in zip(centers, labels):
            axes.text(1.02, value, label, transform=axes.get_yaxis_transform(), horizontalalignment='center')

        axes.yaxis.grid(True, "major")
        axes.xaxis.grid(False, "major")
        axes.legend()
        return axes

    @property
    def observatory(self):
        """
        Retrieves the goes satellite number by parsing the meta dictionary.
        """
        # Various pattern matches for the meta fields.
        pattern_inst = ("GOES 1-{SatelliteNumber:02d} {}")
        pattern_new = ("sci_gxrs-l2-irrad_g{SatelliteNumber:02d}_d{year:4d}{month:2d}{day:2d}_{}.nc")
        pattern_old = ("go{SatelliteNumber:02d}{}{month:2d}{day:2d}.fits")
        pattern_r = ("sci_xrsf-l2-flx1s_g{SatelliteNumber:02d}_d{year:4d}{month:2d}{day:2d}_{}.nc")
        pattern_1m = ("sci_xrsf-l2-avg1m_g{SatelliteNumber:02d}_d{year:4d}{month:2d}{day:2d}_{}.nc")
        pattern_telescop = ("GOES {SatelliteNumber:02d}")
        # The ordering of where we get the metadata from is important.
        # We always want to check ID first as that will most likely have the correct information.
        # The other fields are fallback and sometimes have data in them that is "useless".
        id = (
            self.meta.metas[0].get("id", "").strip()
            or self.meta.metas[0].get("filename_id", "").strip()
            or self.meta.metas[0].get("TELESCOP", "").strip()
            or self.meta.metas[0].get("Instrument", "").strip()
        )
        if isinstance(id, bytes):
            # Needed for h5netcdf < 0.14.0
            id = id.decode('ascii')
        if id is None:
            log.debug("Unable to get a satellite number from 'Instrument', 'TELESCOP' or 'id' ")
            return None
        for pattern in [pattern_inst, pattern_new, pattern_old, pattern_r, pattern_1m, pattern_telescop]:
            parsed = parse(pattern, id)
            if parsed is not None:
                return f"GOES-{parsed['SatelliteNumber']}"
        log.debug('Satellite Number not found in metadata')
        return None

    @classmethod
    def _parse_file(cls, filepath):
        """
        Parses a GOES/XRS FITS file.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        if sunpy.io._file_tools.detect_filetype(filepath) == "hdf5":
            return cls._parse_netcdf(filepath)
        try:
            hdus = sunpy.io._file_tools.read_file(filepath)
        except UnrecognizedFileTypeError:
            raise ValueError(
                f"{Path(filepath).name} is not supported. Only fits and netCDF (nc) can be read.")
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
            # TODO how to extract quality flags from HDU?
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

        # Remove bad values as defined in header comments
        xrsb[xrsb == -99999] = np.nan
        xrsa[xrsa == -99999] = np.nan

        # Fix byte ordering
        newxrsa = xrsa.view(xrsa.dtype.newbyteorder()).byteswap()
        newxrsb = xrsb.view(xrsb.dtype.newbyteorder()).byteswap()

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
        filepath : `str`
            The path of the file to parse
        """
        with h5netcdf.File(filepath, mode="r", **XRSTimeSeries._netcdf_read_kw) as h5nc:
            header = MetaDict(OrderedDict(h5nc.attrs))
            if len(header["id"].strip()) == 0:  # needed to get observatory number if 'id' empty.
                header.update({"filename_id": Path(filepath).name})
            flux_name = h5nc.variables.get("a_flux") or h5nc.variables.get("xrsa_flux")
            if flux_name is None:
                raise ValueError(f"No flux data (either a_flux or xrsa_flux) found in file: {filepath}")
            flux_name_a = flux_name.name
            flux_name_b = flux_name_a.replace("a", "b")
            xrsa = np.asarray(h5nc[flux_name_a])
            xrsb = np.asarray(h5nc[flux_name_b])
            flux_flag_a = h5nc.variables.get("a_flags") or h5nc.variables.get("xrsa_flags") or h5nc.variables.get("xrsa_flag")
            flux_flag_b = h5nc.variables.get("b_flags") or h5nc.variables.get("xrsb_flags") or h5nc.variables.get("xrsb_flag")
            xrsa_quality = np.asarray(h5nc[flux_flag_a.name])
            xrsb_quality = np.asarray(h5nc[flux_flag_b.name])
            start_time_str = h5nc["time"].attrs["units"]
            # h5netcdf < 0.14 return bytes instead of a str
            if isinstance(start_time_str, bytes):
                start_time_str = start_time_str.decode("utf-8")
            start_time_str = start_time_str.lstrip("seconds since").rstrip("UTC").strip()
            times = Time(parse_time(start_time_str).unix + h5nc["time"], format="unix")
            # Checks for primary detector information
            detector_info = False
            if "xrsa_primary_chan" in h5nc:
                detector_info = True
                xrsa_primary_chan = np.asarray(h5nc["xrsa_primary_chan"])
                xrsb_primary_chan = np.asarray(h5nc["xrsb_primary_chan"])
        try:
            times = times.datetime
        except ValueError:
            # We do not make the assumption that the leap second occurs at the end of the file.
            # Therefore, we need to find it:
            # To do so, we convert the times to isot strings, use numpy to find the the leap second string,
            # then use that to workout the index of the leap timestamp.
            idx = np.argwhere((np.char.find(times.isot, ':60.') != -1) == True)  # NOQA: E712
            # We only handle the case there is only 1 leap second in the file.
            # I don't think there every would be a case where it would be more than 1.
            if len(idx) != 1:
                raise ValueError(f"More than one leap second was found in: {Path(filepath).name}")
            warn_user(
                f"There is one leap second timestamp present in: {Path(filepath).name}, "
                "This timestamp has been rounded to `:59.999` to allow its conversion into a Python datetime. "
                f"The leap second timestamp was: {times.isot[idx]}"
            )
            times[idx] = Time(times[idx].isot.tolist()[0][0][:17] + "59.999").unix
            times = times.datetime
        data = DataFrame({"xrsa": xrsa, "xrsb": xrsb, "xrsa_quality": xrsa_quality, "xrsb_quality": xrsb_quality}, index=times)
        units = OrderedDict(
            [
                ("xrsa", u.W/u.m**2),
                ("xrsb", u.W/u.m**2),
                ("xrsa_quality", u.dimensionless_unscaled),
                ("xrsb_quality", u.dimensionless_unscaled),
            ]
        )
        # Adds primary detector info for GOES-R satellites
        if detector_info:
            data["xrsa_primary_chan"] = xrsa_primary_chan
            data["xrsb_primary_chan"] = xrsb_primary_chan
            units.update({"xrsa_primary_chan": u.dimensionless_unscaled,
                          "xrsb_primary_chan": u.dimensionless_unscaled})
        data = data.replace(-9999, np.nan)
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
                if sunpy.io._file_tools.detect_filetype(kwargs["filepath"]) == "hdf5":
                    with h5netcdf.File(kwargs["filepath"], mode="r", **cls._netcdf_read_kw) as f:
                        summary = f.attrs["summary"]
                        if not isinstance(summary, str):
                            # h5netcdf<0.14
                            summary = summary.astype(str)
                        return "XRS" in summary
            except Exception as e:
                log.debug(f'Reading {kwargs["filepath"]} failed with the following exception:\n{e}')
                return False
