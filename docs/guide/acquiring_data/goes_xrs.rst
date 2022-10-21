GOES XRS data
=============

The X-ray Sensor (XRS) on board the GOES series of satellites have provided soft X-ray measurements in two broadband energy ranges 0.5-4 and 1-8 angstrom since 1975.
The GOES 16 and 17 satellites are the latest in line. The flux levels in the GOES 1-8 angstrom channel are used to report flares and determine their size (i.e. their GOES class).

In this guide we are going to look at how you can query and retrieve the GOES XRS data using `~sunpy.net.Fido` and load it into a `~sunpy.timeseries.TimeSeries`.

Some things to note: NOAA have recently re-processed the GOES 13, 14 and 15 XRS science quality data, such that the SWPC scaling factor has been removed.
This means that the fluxes will have a different values, and so will flare peak fluxes from previous 13, 14 and 15 XRS data.
See `here <https://satdat.ngdc.noaa.gov/sem/goes/data/science/xrs/GOES_13-15_XRS_Science-Quality_Data_Readme.pdf>`__ for more details.
The sunpy GOES XRS client for Fido now provides this new re-processed data.
We now also provide the data for GOES 16 and 17.

Another thing to note is that the GOES XRS client `~sunpy.net.Fido` now returns all available GOES data for the specific timerange queried.
For example, there are times when GOES 13, 14 and 15 overlap and such data is available from each satellite.
Similarly there are times when GOES 16 and 17 overlap.

To get started, lets import some packages::

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>>
    >>> from astropy.visualization import time_support
    >>>
    >>> from sunpy import timeseries as ts
    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a

Now lets define our start and end times and query using the `~sunpy.net.Fido`::

    >>> tstart = "2015-06-21 01:00"
    >>> tend = "2015-06-21 23:00"
    >>> result = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"))   # doctest: +REMOTE_DATA
    >>> print(result)   # doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the XRSClient:
    Source: <13: https://umbra.nascom.nasa.gov/goes/fits
    13, 14, 15: https://satdat.ngdc.noaa.gov/sem/goes/data/science/
    16, 17: https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/
    <BLANKLINE>
           Start Time               End Time        Instrument ... Source Provider
    ----------------------- ----------------------- ---------- ... ------ --------
    2015-06-21 00:00:00.000 2015-06-21 23:59:59.999        XRS ...   GOES     NOAA
    2015-06-21 00:00:00.000 2015-06-21 23:59:59.999        XRS ...   GOES     NOAA
    2015-06-21 00:00:00.000 2015-06-21 23:59:59.999        XRS ...   GOES     NOAA

As we can see this now returns three results, one file for GOES 13, one for GOES 14 and one for GOES 15, which can be identified by the ``SatelliteNumber`` column.
However, we probably will only want one of these files for our analysis, so we can query by the `sunpy.net.attrs`: `sunpy.net.dataretriever.attrs.goes.SatelliteNumber` to specify what GOES satellite number we want to use::

    >>> result_goes15 = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"), a.goes.SatelliteNumber(15))  # doctest: +REMOTE_DATA
    >>> print(result_goes15)  # doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the XRSClient:
    Source: <13: https://umbra.nascom.nasa.gov/goes/fits
    13, 14, 15: https://satdat.ngdc.noaa.gov/sem/goes/data/science/
    16, 17: https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/
    <BLANKLINE>
           Start Time               End Time        Instrument ... Source Provider
    ----------------------- ----------------------- ---------- ... ------ --------
    2015-06-21 00:00:00.000 2015-06-21 23:59:59.999        XRS ...   GOES     NOAA


Now we can see that this returns just one file for the GOES 15 data.
Lets now download this data using `~sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch`::

    >>> file_goes15 = Fido.fetch(result_goes15)  # doctest: +REMOTE_DATA

Lets now load this data into a `~sunpy.timeseries.TimeSeries`, and inspect the data by making a plot::

.. plot::

    >>> goes_15 = ts.TimeSeries(file_goes15)  # doctest: +REMOTE_DATA
    >>> fig, ax = plt.subplots()  # doctest: +REMOTE_DATA
    >>> goes_15.plot()  # doctest: +REMOTE_DATA
     <AxesSubplot: ylabel='Watts m$^{-2}$'>

The resulting `~sunpy.timeseries.TimeSeries` can be filtered by GOES quality flags::

    >>> df = goes_15.to_dataframe()  # doctest: +REMOTE_DATA
    >>> df = df[(df["xrsa_quality"] == 0) & (df["xrsb_quality"] == 0)]  # doctest: +REMOTE_DATA
    >>> goes_15 = ts.TimeSeries(df, goes_15.meta, goes_15.units)  # doctest: +REMOTE_DATA

For more information on the flags refer to the `GOES Data Guide <https://satdat.ngdc.noaa.gov/sem/goes/data/science/xrs/GOES_13-15_XRS_Science-Quality_Data_Readme.pdf>`__.

We can also pull out the individual GOES chanels and plot.
The 0.5-4 angstrom channel is known as the "xrsa" channel and the 1-8 angstrom channel is known as the "xrsb" channel::

.. plot::

    >>> fig, ax = plt.subplots()  # doctest: +REMOTE_DATA
    >>> goes_15.plot(columns=["xrsb"])  # doctest: +REMOTE_DATA
    <AxesSubplot: ylabel='W / m2'>

We can also truncate the data for the time of the large flare, and analyze the different channels.
For example, we can plot the derivative which is useful in terms of the Neupert effect when analyzing flares::

.. plot::

    >>> goes_flare = goes_15.truncate("2015-06-21 09:35", "2015-06-21 10:30")   # doctest: +REMOTE_DATA
    >>> time_support()  # doctest: +REMOTE_DATA
    <astropy.visualization.time.time_support.<locals>.MplTimeConverter object at ...>
    >>> fig, ax = plt.subplots()  # doctest: +REMOTE_DATA
    >>> ax.plot(goes_flare.time, np.gradient(goes_flare.quantity("xrsb")))  # doctest: +REMOTE_DATA
    [<matplotlib.lines.Line2D object at ...>]
    >>> ax.set_ylabel("Flux (Wm$^{-2}$$s^{-1}$)")  # doctest: +REMOTE_DATA
    Text(0, 0.5, 'Flux (Wm$^{-2}$$s^{-1}$)')
    >>> fig.autofmt_xdate()  # doctest: +REMOTE_DATA

GOES 16 and 17 data
-------------------
Since March 2020, data prior to GOES 15 (incl) is no longer supported by NOAA and GOES 16 and 17 data is now provided.
See `here <https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/docs/GOES-R_XRS_L2_Data_Users_Guide.pdf>`__ for more details.
GOES 16 and 17 are part of the GOES-R series and provide XRS data at a better time resolution (1s).
sunpy now supports this data also.
GOES 16 has been taking observations from 2017, and GOES 17 since 2018, both of which are now and its now available through sunpy.net.Fido.

Lets query for some recent data over two days::

    >>> results = Fido.search(a.Time("2020-11-20 00:00", "2020-11-21 23:00"), a.Instrument("XRS"))  # doctest: +REMOTE_DATA
    >>> print(results)  # doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    4 Results from the XRSClient:
    Source: <13: https://umbra.nascom.nasa.gov/goes/fits
    13, 14, 15: https://satdat.ngdc.noaa.gov/sem/goes/data/science/
    16, 17: https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/
    <BLANKLINE>
           Start Time               End Time        Instrument ... Source Provider
    ----------------------- ----------------------- ---------- ... ------ --------
    2020-11-20 00:00:00.000 2020-11-20 23:59:59.999        XRS ...   GOES     NOAA
    2020-11-21 00:00:00.000 2020-11-21 23:59:59.999        XRS ...   GOES     NOAA
    2020-11-20 00:00:00.000 2020-11-20 23:59:59.999        XRS ...   GOES     NOAA
    2020-11-21 00:00:00.000 2020-11-21 23:59:59.999        XRS ...   GOES     NOAA

We can see that we are provided with 4 results, two files for GOES 16 and two for GOES 17.
Again we can make the query only specifying one GOES satellite number::

    >>> results_16 = Fido.search(a.Time("2020-11-20 00:00", "2020-11-21 23:00"), a.Instrument("XRS"), a.goes.SatelliteNumber(16))  # doctest: +REMOTE_DATA
    >>> print(results_16)  # doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the XRSClient:
    Source: <13: https://umbra.nascom.nasa.gov/goes/fits
    13, 14, 15: https://satdat.ngdc.noaa.gov/sem/goes/data/science/
    16, 17: https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/
    <BLANKLINE>
           Start Time               End Time        Instrument ... Source Provider
    ----------------------- ----------------------- ---------- ... ------ --------
    2020-11-20 00:00:00.000 2020-11-20 23:59:59.999        XRS ...   GOES     NOAA
    2020-11-21 00:00:00.000 2020-11-21 23:59:59.999        XRS ...   GOES     NOAA

Lets now download this data and load into a `~sunpy.timeseries.TimeSeries`.
We use the `concatenate=True` keyword argument in TimeSeries, as we have two files and want to create one timeseries from them::

.. plot::

    >>> files = Fido.fetch(results_16)  # doctest: +REMOTE_DATA
    >>> goes_16 = ts.TimeSeries(files, concatenate=True)  # doctest: +REMOTE_DATA
    >>> fig, ax = plt.subplots()
    >>> goes_16.plot()  # doctest: +REMOTE_DATA
    <AxesSubplot: ylabel='Watts m$^{-2}$'>
