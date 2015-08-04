# -*- coding: utf-8 -*-
"""Provides programs to process and analyze EVE data."""
from __future__ import absolute_import

import os
import numpy
from datetime import datetime
from collections import OrderedDict

import matplotlib.pyplot as plt
from pandas.io.parsers import read_csv
from os.path import basename

from sunpy.lightcurve import LightCurve

__all__ = ['EVELightCurve']


class EVELightCurve(LightCurve):
    """
    SDO EVE LightCurve for level 0CS data.

    The Extreme Ultraviolet Variability Experiment (EVE) is an instrument on board
    the Solar Dynamics Observatory (SDO). The EVE instrument is designed to
    measure the solar extreme ultraviolet (EUV) irradiance. The EUV radiation
    includes the 0.1-105 nm range, which provides the majority of the energy
    for heating Earth’s thermosphere and creating Earth’s ionosphere (charged plasma).

    EVE includes several irradiance instruments: The Multiple EUV Grating
    Spectrographs (MEGS)-A is a grazing- incidence spectrograph that measures
    the solar EUV irradiance in the 5 to 37 nm range with 0.1-nm resolution,
    and the MEGS-B is a normal-incidence, dual-pass spectrograph that measures
    the solar EUV irradiance in the 35 to 105 nm range with 0.1-nm resolution.

    Level 0CS data is primarily used for space weather. It is provided near
    real-time and is crudely calibrated 1-minute averaged broadband irradiances
    from ESP and MEGS-P broadband.

    Data is available starting on 2010/03/01.

    Examples
    --------
    >>> import sunpy.lightcurve
    >>> import sunpy.data.test

    >>> eve = sunpy.lightcurve.EVELightCurve.create()
    >>> eve = sunpy.lightcurve.EVELightCurve.create('2012/06/20')
    >>> eve = sunpy.lightcurve.EVELightCurve.create(sunpy.data.test.EVE_AVERAGES_CSV)
    >>> eve = sunpy.lightcurve.EVELightCurve.create("http://lasp.colorado.edu/eve/data_access/quicklook/quicklook_data/L0CS/LATEST_EVE_L0CS_DIODES_1m.txt")
    >>> eve.peek(subplots=True)    # doctest: +SKIP

    References
    ----------
    * `SDO Mission Homepage <http://sdo.gsfc.nasa.gov>`_
    * `EVE Homepage <http://lasp.colorado.edu/home/eve/>`_
    * `Level 0CS Definition <http://lasp.colorado.edu/home/eve/data/>`_
    * `EVE Data Acess <http://lasp.colorado.edu/home/eve/data/data-access/>`_
    * `Instrument Paper <http://link.springer.com/article/10.1007%2Fs11207-009-9487-6>`_
    """

    def peek(self, column=None, **kwargs):
        """Plots the light curve in a new figure. An example is shown below.

        .. plot::

            import sunpy
            eve = sunpy.lightcurve.EVELightCurve.create('2012/06/20')
            eve.peek(subplots=True)

        Parameters
        ----------
        column : str
            The column to display. If None displays all.

        **kwargs : dict
            Any additional plot arguments that should be used
            when plotting.

        Returns
        -------
        fig : `~matplotlib.Figure`
            A plot figure.
        """
        figure = plt.figure()
        # Choose title if none was specified
        if "title" not in kwargs and column is None:
            if len(self.data.columns) > 1:
                kwargs['title'] = 'EVE (1 minute data)'
            else:
                if self._filename is not None:
                    base = self._filename.replace('_', ' ')
                    kwargs['title'] = os.path.splitext(base)[0]
                else:
                    kwargs['title'] = 'EVE Averages'

        if column is None:
            self.plot(**kwargs)
        else:
            data = self.data[column]
            if "title" not in kwargs:
                kwargs['title'] = 'EVE ' + column.replace('_', ' ')
            data.plot(**kwargs)
        figure.show()
        return figure

    @staticmethod
    def _get_default_uri():
        """Loads latest level 0CS if no other data is specified"""
        return "http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/LATEST_EVE_L0CS_DIODES_1m.txt"

    @staticmethod
    def _get_url_for_date(date):
        """Returns a URL to the EVE data for the specified date

            @NOTE: currently only supports downloading level 0 data
            .TODO: No data available prior to 2010/03/01!
        """
        base_url = 'http://lasp.colorado.edu/eve/data_access/evewebdata/quicklook/L0CS/SpWx/'
        return base_url + date.strftime('%Y/%Y%m%d') + '_EVE_L0CS_DIODES_1m.txt'

    @classmethod
    def _parse_csv(cls, filepath):
        """Parses an EVE CSV file."""
        cls._filename = basename(filepath)
        with open(filepath, 'rb') as fp:
            # Determine type of EVE CSV file and parse
            line1 = fp.readline()
            fp.seek(0)

            if line1.startswith("Date"):
                return cls._parse_average_csv(fp)
            elif line1.startswith(";"):
                return cls._parse_level_0cs(fp)

    @staticmethod
    def _parse_average_csv(fp):
        """Parses an EVE Averages file."""
        return "", read_csv(fp, sep=",", index_col=0, parse_dates=True)

    @staticmethod
    def _parse_level_0cs(fp):
        """Parses and EVE Level 0CS file."""
        is_missing_data = False      #boolean to check for missing data
        missing_data_val = numpy.nan
        header = []
        fields = []
        line = fp.readline()
        # Read header at top of file
        while line.startswith(";"):
            header.append(line)
            if '; Missing data:' in line :
                is_missing_data = True
                missing_data_val = line.split(':')[1].strip()

            line = fp.readline()

        meta = OrderedDict()
        for hline in header :
            if hline == '; Format:\n' or hline == '; Column descriptions:\n':
                continue
            elif ('Created' in hline) or ('Source' in hline):
                meta[hline.split(':',1)[0].replace(';',' ').strip()] = hline.split(':',1)[1].strip()
            elif ':' in hline :
                meta[hline.split(':')[0].replace(';',' ').strip()] = hline.split(':')[1].strip()

        fieldnames_start = False
        for hline in header:
            if hline.startswith("; Format:"):
                fieldnames_start = False
            if fieldnames_start:
                fields.append(hline.split(":")[0].replace(';', ' ').strip())
            if hline.startswith("; Column descriptions:"):
                fieldnames_start = True

        # Next line is YYYY DOY MM DD
        date_parts = line.split(" ")

        year = int(date_parts[0])
        month = int(date_parts[2])
        day = int(date_parts[3])
        #last_pos = fp.tell()
        #line = fp.readline()
        #el = line.split()
        #len

        # function to parse date column (HHMM)
        parser = lambda x: datetime(year, month, day, int(x[0:2]), int(x[2:4]))

        data = read_csv(fp, sep="\s*", names=fields, index_col=0, date_parser=parser, header = None)
        if is_missing_data :   #If missing data specified in header
            data[data == float(missing_data_val)] = numpy.nan

        # data.columns = fields
        return meta, data
