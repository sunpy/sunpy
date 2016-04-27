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

from sunpy import config
TIME_FORMAT = config.get("general", "time_format")

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

    def plot(self, title=True, axes=None, plot_type=None, **plot_args):
        """Plots EVE light curve is the usual manner"""
        if axes is None:
            axes = plt.gca()

        if plot_type == None:
            plot_type = self._get_plot_types()[0]

        if plot_type == self._get_plot_types()[0]:      # goes proxy
            self.data['XRS-B proxy'].plot(ax=axes, label='1.0--8.0 $\AA$', color='red', lw=2, **plot_args)
            self.data['XRS-A proxy'].plot(ax=axes, label='0.5--4.0 $\AA$', color='blue', lw=2, **plot_args)
            axes.set_yscale("log")
            axes.set_ylim(1e-9, 1e-2)
            axes.set_ylabel(self.meta.get('UNIT')[0])
            ax2 = axes.twinx()
            ax2.set_yscale("log")
            ax2.set_ylim(1e-9, 1e-2)
            ax2.set_yticks((1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
            ax2.set_yticklabels((' ', 'A', 'B', 'C', 'M', 'X', ' '))
        elif plot_type == self._get_plot_types()[1]:  # esp quad
            self.data['q0ESP'].plot(ax=axes, label='ESP 0', **plot_args)
            self.data['q1ESP'].plot(ax=axes, label='ESP 1', **plot_args)
            self.data['q2ESP'].plot(ax=axes, label='ESP 2', **plot_args)
            self.data['q3ESP'].plot(ax=axes, label='ESP 3', **plot_args)
        elif plot_type == self._get_plot_types()[2]:  # position
            self.data['CMLat'].plot(ax=axes, label='Latitude', **plot_args)
            self.data['CMLon'].plot(ax=axes, label='Longitude', **plot_args)
        elif plot_type == self._get_plot_types()[3]:  # sem
            self.data['SEM proxy'].plot(ax=axes, label='SEM Proxy', **plot_args)
        elif plot_type == self._get_plot_types()[4]:  # esp
            self.data['17.1ESP'].plot(ax=axes, label='17.1', **plot_args)
            self.data['25.7ESP'].plot(ax=axes, label='25.7', **plot_args)
            self.data['30.4ESP'].plot(ax=axes, label='30.4', **plot_args)
            self.data['36.6ESP'].plot(ax=axes, label='36.6', **plot_args)
        elif plot_type == self._get_plot_types()[5]:  # dark
            self.data['darkESP'].plot(ax=axes, label='ESP', **plot_args)
            self.data['darkMEGS-P'].plot(ax=axes, label='MEGS-P', **plot_args)
        else:
            raise ValueError('Not a recognized plot type.')

        axes.set_xlabel('Start time: ' + self.data.index[0].strftime(TIME_FORMAT))
        axes.set_title(title)
        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        plt.gcf().autofmt_xdate()
        return axes

    @classmethod
    def _get_plot_types(cls):
        return ['goes proxy', 'esp quad', 'position', 'sem', 'esp', 'dark']

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

            if line1.startswith("Date".encode('ascii')):
                return cls._parse_average_csv(fp)
            elif line1.startswith(";".encode('ascii')):
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
        while line.startswith(";".encode('ascii')):
            header.append(line)
            if '; Missing data:'.encode('ascii') in line :
                is_missing_data = True
                missing_data_val = line.split(':'.encode('ascii'))[1].strip()

            line = fp.readline()

        meta = OrderedDict()
        for hline in header :
            if hline == '; Format:\n'.encode('ascii') or hline == '; Column descriptions:\n'.encode('ascii'):
                continue
            elif ('Created'.encode('ascii') in hline) or ('Source'.encode('ascii') in hline):
                meta[hline.split(':'.encode('ascii'),
                                 1)[0].replace(';'.encode('ascii'),
                                               ' '.encode('ascii')).strip()] = hline.split(':'.encode('ascii'), 1)[1].strip()
            elif ':'.encode('ascii') in hline :
                meta[hline.split(':'.encode('ascii'))[0].replace(';'.encode('ascii'), ' '.encode('ascii')).strip()] = hline.split(':'.encode('ascii'))[1].strip()

        fieldnames_start = False
        for hline in header:
            if hline.startswith("; Format:".encode('ascii')):
                fieldnames_start = False
            if fieldnames_start:
                fields.append(hline.split(":".encode('ascii'))[0].replace(';'.encode('ascii'), ' '.encode('ascii')).strip())
            if hline.startswith("; Column descriptions:".encode('ascii')):
                fieldnames_start = True

        # Next line is YYYY DOY MM DD
        date_parts = line.split(" ".encode('ascii'))

        year = int(date_parts[0])
        month = int(date_parts[2])
        day = int(date_parts[3])
        #last_pos = fp.tell()
        #line = fp.readline()
        #el = line.split()
        #len

        # function to parse date column (HHMM)
        parser = lambda x: datetime(year, month, day, int(x[0:2]), int(x[2:4]))

        data = read_csv(fp, sep="\s*".encode('ascii'), names=fields, index_col=0, date_parser=parser, header=None, engine='python')
        if is_missing_data :   #If missing data specified in header
            data[data == float(missing_data_val)] = numpy.nan

        # data.columns = fields
        units = [None] * len(data.columns)
        units[0] = units[1] = 'Watts m^-2'
        meta.update({'unit': units})
        meta.update({'instrume': 'SDO/EVE'})
        meta.update({'obsrvtry': 'SDO'})
        meta.update({'telescope': 'EVE'})
        meta.update({'wavelnth': ''})
        meta.update({'waveunit': ''})
        return meta, data
