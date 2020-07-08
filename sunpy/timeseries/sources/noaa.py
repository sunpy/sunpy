"""
This module provies NOAA Solar Cycle `~sunpy.timeseries.TimeSeries` source.
"""
from pathlib import Path
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv

import astropy.units as u
from astropy.time import Time

from sunpy.time import parse_time
from sunpy.timeseries.timeseriesbase import GenericTimeSeries
from sunpy.util.decorators import deprecated
from sunpy.util.metadata import MetaDict
from sunpy.visualization import peek_show

__all__ = ['NOAAIndicesTimeSeries', 'NOAAPredictIndicesTimeSeries']


class NOAAIndicesTimeSeries(GenericTimeSeries):
    """
    NOAA Solar Cycle monthly indices.

    Solar activity is measured by a number of different values.
    The NOAA Solar Weather Prediction Center (SWPC) publishes the following indices.
    All of these indices are also provided as a 13-month running smoothed value.

    * The SWO sunspot number is issued by the NOAA Space Weather Prediction Center (SWPC)
    * The RI sunspot number is the official International Sunspot Number and is
      issued by the `Solar Influence Data Analysis Center (SDIC) <http://sidc.oma.be>`_ in Brussels, Belgium.
    * The ratio between the SWO and RI indices.
    * Radio flux at 10.7 cm is produced by `Penticon/Ottawa <https://www.ngdc.noaa.gov/stp/solar/flux.html>`_ and the units are in sfu.
    * The Ap Geomagnetic Index is produced by the United States Air Force (USAF).

    .. note::
        See the gallery example :ref:`sphx_glr_generated_gallery_plotting_solar_cycle_example.py`
        for how to use `~sunpy.net.Fido` to retrieve the data file.

    Examples
    --------
    >>> import sunpy.timeseries
    >>> noaa_url = "https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json"
    >>> noaa = sunpy.timeseries.TimeSeries(noaa_url, source='NOAAIndices')  # doctest: +REMOTE_DATA
    >>> noaa.peek()  # doctest: +SKIP

    References
    ----------
    * `Solar and Geomagnetic Indices Data Archive <https://www.swpc.noaa.gov/products/3-day-geomagnetic-forecast>`_
    * `Recent solar indices <https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json>`_
    * `Indices Descriptions <ftp://ftp.swpc.noaa.gov/pub/weekly/README3>`_
    * `NOAA plots of Solar Cycle Progression <https://www.swpc.noaa.gov/products/solar-cycle-progression>`_
    * `NOAA Product List <https://www.swpc.noaa.gov/products-and-data>`_
    """
    # Class attribute used to specify the source class of the TimeSeries.
    _source = 'noaaindices'

    @peek_show
    def peek(self, type='sunspot SWO', **kwargs):
        """
        Plots NOAA Indices as a function of time. An example is shown below.

        .. plot::

            import sunpy.timeseries
            noaa_url = "https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json"
            noaa = sunpy.timeseries.TimeSeries(noaa_url, source='NOAAIndices')
            noaa.peek()

        Parameters
        ----------
        type : {'sunspot SWO', 'sunspot RI', 'sunspot compare', 'radio', 'geo'}, optional
            The type of plot required. Defaults to "sunspot SWO".
        **kwargs : `dict`
            Additional plot keyword arguments that are handed to `axes.plot` functions
        """
        # Check we have a timeseries valid for plotting
        self._validate_data_for_plotting()

        figure = plt.figure()
        axes = plt.gca()
        dataframe = self.to_dataframe()
        if "sunspot SWO" in dataframe.columns.values:
            if type == 'sunspot SWO':
                axes = dataframe['sunspot SWO'].plot(**kwargs)
                dataframe['sunspot SWO smooth'].plot(**kwargs)
                axes.set_ylabel('Sunspot Number')
            elif type == 'sunspot RI':
                axes = dataframe['sunspot RI'].plot(**kwargs)
                dataframe['sunspot RI smooth'].plot(**kwargs)
                axes.set_ylabel('Sunspot Number')
            elif type == 'sunspot compare':
                axes = dataframe['sunspot RI'].plot(**kwargs)
                dataframe['sunspot SWO'].plot(**kwargs)
                axes.set_ylabel('Sunspot Number')
            elif type == 'radio':
                axes = dataframe['radio flux'].plot(**kwargs)
                dataframe['radio flux smooth'].plot(**kwargs)
                axes.set_ylabel('Radio Flux [sfu]')
            elif type == 'geo':
                axes = dataframe['geomagnetic ap'].plot(**kwargs)
                dataframe['geomagnetic ap smooth'].plot(**kwargs)
                axes.set_ylabel('Geomagnetic AP Index')
            else:
                raise ValueError(f'Got unknown plot type "{type}"')
        else:
            dataframe['mean monthly S.I.D.C. sunspot number'].plot(**kwargs)
            dataframe['smoothed S.I.D.C. sunspot number'].plot(**kwargs)
            dataframe['mean monthly SWPC/SWO sunspot number'].plot(**kwargs)
            dataframe['smoothed SWPC/SWO sunspot number'].plot(**kwargs)
            dataframe['mean monthly 10.7cm radio flux'].plot(**kwargs)
            dataframe['smoothed 10.7cm radio flux'].plot(**kwargs)

        axes.set_ylim(0)
        axes.set_title('Solar Cycle Progression')

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.legend()

        return figure

    @classmethod
    def _parse_file(cls, filepath):
        suffix = Path(filepath).suffix
        if suffix == '.json':
            return cls._parse_json_file(filepath)
        elif suffix == ".txt":
            return cls._parse_txt_file(filepath)
        else:
            raise ValueError(f"{Path(filepath).name} does not have a suffix of '.txt' or '.json'")

    @staticmethod
    @deprecated("2.1", "NOAA data products have moved to a new JSON file format.")
    def _parse_txt_file(filepath):
        """
        Parses an NOAA indices text file.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        header = []
        with open(filepath, 'r') as fp:
            line = fp.readline()
            # Read header at top of file
            while line.startswith((":", "#")):
                header += line
                line = fp.readline()
            fields = ('yyyy', 'mm', 'sunspot SWO', 'sunspot RI', 'sunspot ratio', 'sunspot SWO smooth',
                      'sunspot RI smooth', 'radio flux', 'radio flux smooth', 'geomagnetic ap', 'geomagnetic smooth')
            data = read_csv(fp, delim_whitespace=True, names=fields,
                            comment='#', dtype={'yyyy': np.str, 'mm': np.str})
            data = data.dropna(how='any')
            timeindex = Time.strptime([x + y for x, y in zip(data['yyyy'], data['mm'])], '%Y%m')
            timeindex.precision = 9
            data['time'] = timeindex.isot.astype('datetime64')
            data = data.set_index('time')
            data = data.drop('mm', 1)
            data = data.drop('yyyy', 1)

            # Add the units data
            units = OrderedDict([('sunspot SWO', u.dimensionless_unscaled),
                                 ('sunspot RI', u.dimensionless_unscaled),
                                 ('sunspot ratio', u.dimensionless_unscaled),
                                 ('sunspot SWO smooth', u.dimensionless_unscaled),
                                 ('sunspot RI smooth', u.dimensionless_unscaled),
                                 ('radio flux', u.W/u.m**2),
                                 ('radio flux smooth', u.W/u.m**2),
                                 ('geomagnetic ap', u.dimensionless_unscaled),
                                 ('geomagnetic smooth', u.dimensionless_unscaled)])
            # TODO: check units
            # TODO: fix header/meta, it's returning rubbish.
            return data, MetaDict({'comments': header}), units

    @staticmethod
    def _parse_json_file(filepath):
        """
        Parses an NOAA indices JSON file.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        with open(filepath) as fp:
            fp.seek(0)
            data = pd.read_json(fp.read())

        rename = {'ssn': 'mean monthly S.I.D.C. sunspot number',
                  'smoothed_ssn': 'smoothed S.I.D.C. sunspot number',
                  'observed_swpc_ssn': 'mean monthly SWPC/SWO sunspot number',
                  'smoothed_swpc_ssn': 'smoothed SWPC/SWO sunspot number',
                  'f10.7': 'mean monthly 10.7cm radio flux',
                  'smoothed_f10.7': 'smoothed 10.7cm radio flux'}
        data = data.rename(columns=rename)
        # Convoluted time index handling
        data = data.set_index('time-tag')
        data.index = pd.DatetimeIndex(data.index.values)
        data.index = pd.DatetimeIndex(parse_time(
            [x for x in data.index.values]).isot.astype('datetime64'))

        # Add the units data, reported in radio flux values (sfu) originally.
        units = OrderedDict([('mean monthly S.I.D.C. sunspot number', u.dimensionless_unscaled),
                             ('smoothed S.I.D.C. sunspot number', u.dimensionless_unscaled),
                             ('mean monthly SWPC/SWO sunspot number', u.dimensionless_unscaled),
                             ('smoothed SWPC/SWO sunspot number', u.dimensionless_unscaled),
                             ('mean monthly 10.7cm radio flux', 1e-22*u.W/(u.m**2*u.Hertz)),
                             ('smoothed 10.7cm radio flux', 1e-22*u.W/(u.m**2*u.Hertz))])
        return data, MetaDict({'comments': ""}), units

    @classmethod
    def is_datasource_for(cls, **kwargs):
        """
        Determines if header corresponds to an NOAA indices timeseries.
        """
        if kwargs.get('source', ''):
            return kwargs.get('source', '').lower().startswith(cls._source)


class NOAAPredictIndicesTimeSeries(GenericTimeSeries):
    """
    NOAA Solar Cycle Predicted Progression.

    The predictions are updated monthly and are produced by ISES. Observed
    values are initially the preliminary values which are replaced with the
    final values as they become available.

    The following predicted values are available.

    * The predicted RI sunspot number is the official International Sunspot
      Number and is issued by the `Solar Influence Data Analysis Center (SDIC) <http://sidc.oma.be>`_ in Brussels, Belgium.
    * The predicted radio flux at 10.7 cm is produced by
      `Penticon/Ottawa <https://www.ngdc.noaa.gov/stp/solar/flux.html>`_ and the units are in sfu.

    .. note::
        See the gallery example :ref:`sphx_glr_generated_gallery_plotting_solar_cycle_example.py`
        for how to use `~sunpy.net.Fido` to retrieve the data file.

    Examples
    --------
    >>> import sunpy.timeseries
    >>> noaa_url = 'https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json'  # doctest: +REMOTE_DATA
    # doctest: +REMOTE_DATA
    >>> noaa = sunpy.timeseries.TimeSeries(noaa_url, source='NOAAPredictIndices')
    >>> noaa.peek()  # doctest: +SKIP

    References
    ----------
    * `Solar and Geomagnetic Indices Data Archive <https://www.swpc.noaa.gov/products/3-day-geomagnetic-forecast>`_
    * `Predicted solar indices <https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json>`_
    * `NOAA plots of Solar Cycle Progression <https://www.swpc.noaa.gov/products/solar-cycle-progression>`_
    * `NOAA Product List <https://www.swpc.noaa.gov/products-and-data>`_
    """

    # Class attribute used to specify the source class of the TimeSeries.
    _source = 'noaapredictindices'

    @ peek_show
    def peek(self, **plot_args):
        """
        Plots predicted NOAA Indices as a function of time. An example is shown
        below.

        .. plot::

            import sunpy.timeseries
            noaa_url = 'https://services.swpc.noaa.gov/json/solar-cycle/predicted-solar-cycle.json'
            noaa = sunpy.timeseries.TimeSeries(noaa_url, source='NOAAPredictIndices')
            noaa.peek()

        Parameters
        ----------
        **plot_args : `dict`
            Additional plot keyword arguments that are handed to `axes.plot` functions
        """
        # Check we have a timeseries valid for plotting
        self._validate_data_for_plotting()
        figure = plt.figure()
        axes = plt.gca()
        # Old dataframe vs New dataframe
        dataframe = self.to_dataframe()
        if "sunspot" in dataframe.columns.values:
            axes = dataframe['sunspot'].plot(color='b', **plot_args)
            dataframe['sunspot low'].plot(linestyle='--', color='b', **plot_args)
            dataframe['sunspot high'].plot(linestyle='--', color='b', **plot_args)
        else:
            axes = dataframe['predicted sunspot number'].plot(color='b', **plot_args)
            dataframe['predicted sunspot number high range'].plot(
                linestyle='--', color='b', **plot_args)
            dataframe['predicted sunspot number low range'].plot(
                linestyle='--', color='b', **plot_args)

        axes.set_ylim(0)
        axes.set_title('Solar Cycle Sunspot Number Prediction')
        axes.set_ylabel('Sunspot Number')

        axes.yaxis.grid(True, 'major')
        axes.xaxis.grid(True, 'major')
        axes.legend()

        return figure

    @ classmethod
    def _parse_file(cls, filepath):
        suffix = Path(filepath).suffix
        if suffix == '.json':
            return cls._parse_json_file(filepath)
        elif suffix == ".txt":
            return cls._parse_txt_file(filepath)
        else:
            raise ValueError(f"{Path(filepath).name} does not have a suffix of '.txt' or '.json'")

    @ classmethod
    def is_datasource_for(cls, **kwargs):
        """
        Determines if header corresponds to an NOAA predict indices
        `~sunpy.timeseries.TimeSeries`.
        """
        if kwargs.get('source', ''):
            return kwargs.get('source', '').lower().startswith(cls._source)

    @ staticmethod
    @ deprecated("2.1", "NOAA data products have moved to a new JSON file format.")
    def _parse_txt_file(filepath):
        """
        Parses an NOAA Predict indices text file.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        header = ''
        with open(filepath, 'r') as fp:
            line = fp.readline()
            # Read header at top of file
            while line.startswith((":", "#")):
                header += line
                line = fp.readline()
            fields = ('yyyy', 'mm', 'sunspot', 'sunspot low', 'sunspot high',
                      'radio flux', 'radio flux low', 'radio flux high')
            data = read_csv(filepath, delim_whitespace=True, names=fields,
                            comment='#', skiprows=2, dtype={'yyyy': np.str, 'mm': np.str})
            data = data.dropna(how='any')

            timeindex = Time.strptime([x + y for x, y in zip(data['yyyy'], data['mm'])], '%Y%m')
            timeindex.precision = 9
            data['time'] = timeindex.isot.astype('datetime64')

            data = data.set_index('time')
            data = data.drop('mm', 1)
            data = data.drop('yyyy', 1)

            # Add the units data
            units = OrderedDict([('sunspot', u.dimensionless_unscaled),
                                 ('sunspot low', u.dimensionless_unscaled),
                                 ('sunspot high', u.dimensionless_unscaled),
                                 ('radio flux', u.W/u.m**2),
                                 ('radio flux low', u.W/u.m**2),
                                 ('radio flux high', u.W/u.m**2)])
            # Todo: check units used.
            return data, MetaDict({'comments': header}), units

    @ staticmethod
    def _parse_json_file(filepath):
        """
        Parses an NOAA Predict indices JSON file.

        Parameters
        ----------
        filepath : `str`
            The path to the file you want to parse.
        """
        with open(filepath) as fp:
            fp.seek(0)
            data = pd.read_json(fp.read())
        rename = {'predicted_ssn': 'predicted sunspot number',
                  'high_ssn': 'predicted sunspot number high range',
                  'low_ssn': 'predicted sunspot number low range',
                  'predicted_f10.7': 'predicted f10.7cm',
                  'high_f10.7': 'predicted f10.7cm high range',
                  'low_f10.7': 'predicted f10.7cm low range'}
        data = data.rename(columns=rename)
        # Convoluted time index handling
        data = data.set_index('time-tag')
        data.index = pd.DatetimeIndex(data.index.values)
        data.index = pd.DatetimeIndex(parse_time(
            [x for x in data.index.values]).isot.astype('datetime64'))
        # Add the units data, reported in radio flux values (sfu) originally.
        units = OrderedDict([('predicted sunspot number', u.dimensionless_unscaled),
                             ('predicted sunspot number high range', u.dimensionless_unscaled),
                             ('predicted sunspot number low range', u.dimensionless_unscaled),
                             ('predicted f10.7cm', 1e-22*u.W/(u.m**2*u.Hertz)),
                             ('predicted f10.7cm high range', 1e-22*u.W/(u.m**2*u.Hertz)),
                             ('predicted f10.7cm low range', 1e-22*u.W/(u.m**2*u.Hertz))])
        return data, MetaDict({'comments': ""}), units
