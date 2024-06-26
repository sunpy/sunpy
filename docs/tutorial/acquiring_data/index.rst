.. _sunpy-tutorial-acquiring-data-index:

**************
Acquiring Data
**************

This section of the tutorial introduces ways to obtain different kind of solar data from different places.
The main tutorial below focuses on ``Fido``, which is a generic search interface that sunpy provides.
There are some data-provider specific tutorials that you might want to visit after going through this tutorial:

.. toctree::
    :maxdepth: 1

    jsoc
    hek

This guide outlines how to search for and download data using the `~sunpy.net.Fido` interface for search and download.
`~sunpy.net.Fido` is a unified interface for searching and fetching solar physics data irrespective of the underlying client or web service through which the data is obtained.
It therefore supplies a single, easy, and consistent way to obtain most forms of solar physics data.

The `~sunpy.net.Fido` object is in `sunpy.net`.
All the examples in this guide use ``Fido``, so lets start by importing it:

.. code-block:: python

    >>> from sunpy.net import Fido, attrs as a

.. contents::
    :depth: 3

Fido supports a number of different remote data sources. To see a list the Fido object can be printed:

.. code-block:: python

    >>> print(Fido)
    sunpy.net.Fido
    <BLANKLINE>
    Fido is a unified data search and retrieval tool.
    <BLANKLINE>
    It provides simultaneous access to a variety of online data sources, some
    cover multiple instruments and data products like the Virtual Solar
    Observatory and some are specific to a single source.
    <BLANKLINE>
    For details of using `~sunpy.net.Fido` see :ref:`sunpy-tutorial-acquiring-data-index`.
    <BLANKLINE>
    <BLANKLINE>
          Client                                                            Description
    ----------------- -----------------------------------------------------------------------------------------------------------------------
    CDAWEBClient      Provides access to query and download from the Coordinated Data Analysis Web (CDAWeb).
    ADAPTClient       Provides access to the ADvanced Adaptive Prediction Technique (ADAPT) products of the National Solar Observatory (NSO).
    EVEClient         Provides access to Level 0CS Extreme ultraviolet Variability Experiment (EVE) data.
    GBMClient         Provides access to data from the Gamma-Ray Burst Monitor (GBM) instrument on board the Fermi satellite.
    XRSClient         Provides access to several GOES XRS files archive.
    SUVIClient        Provides access to data from the GOES Solar Ultraviolet Imager (SUVI).
    GONGClient        Provides access to the Magnetogram products of NSO-GONG synoptic Maps.
    LYRAClient        Provides access to the LYRA/Proba2 data archive.
    NOAAIndicesClient Provides access to the NOAA solar cycle indices.
    NOAAPredictClient Provides access to the NOAA SWPC predicted sunspot Number and 10.7 cm radio flux values.
    SRSClient         Provides access to the NOAA SWPC solar region summary data.
    NoRHClient        Provides access to the Nobeyama RadioHeliograph (NoRH) averaged correlation time series data.
    RHESSIClient      Provides access to the RHESSI observing summary time series data.
    HEKClient         Provides access to the Heliophysics Event Knowledgebase (HEK).
    HECClient         Provides access to the HELIO webservices.
    JSOCClient        Provides access to the JSOC Data Export service.
    VSOClient         Provides access to query and download from Virtual Solar Observatory (VSO).

Searching for Data
******************

To search for data with `~sunpy.net.Fido`, you need to specify attributes to search with.
Examples of generic search attributes that work across many different data sources are:

- `a.Time <sunpy.net.attrs.Time>`
- `a.Instrument <sunpy.net.attrs.Instrument>`
- `a.Wavelength <sunpy.net.attrs.Wavelength>`

Some other attributes are client specific, and are found under client specific submodules, e.g. `attrs.vso <sunpy.net.vso.attrs>` and `attrs.jsoc <sunpy.net.jsoc.attrs>`.
The full list of attributes can be found in the `attrs submodule reference <sunpy.net.attrs>`.

Some search attributes need one or more values specifying, for example ``Time`` needs at least a start and an end date to specify a time range:

.. code-block:: python

    >>> a.Time('2012/3/4', '2012/3/6')
    <sunpy.net.attrs.Time(2012-03-04 00:00:00.000, 2012-03-06 00:00:00.000)>

For attributes that can take a range of different values, printing the attribute lists the values sunpy knows about.
These values are updated with every release of sunpy, so may not be always up to date!
As an example:

.. code-block:: python

    >>> print(a.Instrument)
    sunpy.net.attrs.Instrument
    <BLANKLINE>
    Specifies the Instrument name for the search.
    <BLANKLINE>
           Attribute Name          Client   ...                                   Description
    --------------------------- ----------- ... --------------------------------------------------------------------------------
    adapt                       ADAPT       ... ADvanced Adaptive Prediction Technique.
    aia                         VSO         ... Atmospheric Imaging Assembly
    bcs                         VSO         ... Bragg Crystal Spectrometer
    be_continuum                VSO         ... INAF-OACT Barra Equatoriale Continuum Instrument
    be_halpha                   VSO         ... INAF-OACT Barra Equatoriale Hα Instrument
    bigbear                     VSO         ... Big Bear Solar Observatory, California TON and GONG+ sites
    ...
    xrs                         XRS         ... GOES X-ray Sensor
    xrt                         VSO         ... X-Ray Telescope

This is a full list of known values, a description, and which clients support those values (if you want to search using a specific data source).
Printing attributes like this is supported for most attributes, including  client specific ones.
These attributes also support tab-completion to auto to auto-fill the attribute name, for example typing ``a.jsoc.aia_f<TAB>`` in a jupyter notebook will show you the available attributes that start with "aia_f".

To search for data use the ``Fido.search`` method:

.. code-block:: python

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.lyra, a.Level.two) # doctest: +REMOTE_DATA

this returns an `~sunpy.net.fido_factory.UnifiedResponse` object containing all the search results that match the search attributes.
This does not download the files; we'll learn how to do that later in :ref:`sunpy-tutorial-acquiring-data-downloading-data`.

To see a summary of the results print the result variable that came back from the previous search:

.. code-block:: python

    >>> print(result)  # doctest: +REMOTE_DATA
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the LYRAClient:
    Source: http://proba2.oma.be/lyra/data/bsd
    <BLANKLINE>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999       LYRA ...      ESA     2
    2012-03-05 00:00:00.000 2012-03-05 23:59:59.999       LYRA ...      ESA     2
    2012-03-06 00:00:00.000 2012-03-06 23:59:59.999       LYRA ...      ESA     2
    <BLANKLINE>
    <BLANKLINE>

Queries can be made more flexible or specific by adding more attrs objects to the `~sunpy.net.Fido` search.
As an example, specific passbands can be searched for by supplying a `~astropy.units.Quantity` to the `a.Wavelength <sunpy.net.attrs.Wavelength>` attribute:

.. code-block:: python

    >>> import astropy.units as u

    >>> Fido.search(a.Time('2012/3/4', '2012/3/4'), a.Instrument.aia,
    ...             a.Wavelength(171*u.angstrom))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 67.789 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2012-03-04 00:00:00.000 2012-03-04 00:00:01.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    <BLANKLINE>

Data of a given cadence can also be specified using the `a.Sample <sunpy.net.attrs.Sample>` attribute:

.. code-block:: python

    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.aia,
    ...             a.Wavelength(171*u.angstrom), a.Sample(120*u.minute))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    25 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 1.695 Gbyte
    <BLANKLINE>
           Start Time               End Time        Source Instrument   Wavelength   Provider  Physobs  Wavetype Extent Width Extent Length Extent Type   Size
                                                                         Angstrom                                                                        Mibyte
    ----------------------- ----------------------- ------ ---------- -------------- -------- --------- -------- ------------ ------------- ----------- --------
    2012-03-04 00:00:00.000 2012-03-04 00:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-04 02:00:00.000 2012-03-04 02:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-04 04:00:00.000 2012-03-04 04:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-04 06:00:00.000 2012-03-04 06:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-04 08:00:00.000 2012-03-04 08:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-04 10:00:00.000 2012-03-04 10:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-04 12:00:00.000 2012-03-04 12:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-04 14:00:00.000 2012-03-04 14:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-04 16:00:00.000 2012-03-04 16:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
                        ...                     ...    ...        ...            ...      ...       ...      ...          ...           ...         ...      ...
    2012-03-05 06:00:00.000 2012-03-05 06:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-05 08:00:00.000 2012-03-05 08:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-05 10:00:00.000 2012-03-05 10:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-05 12:00:00.000 2012-03-05 12:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-05 14:00:00.000 2012-03-05 14:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-05 16:00:00.000 2012-03-05 16:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-05 18:00:00.000 2012-03-05 18:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-05 20:00:00.000 2012-03-05 20:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-05 22:00:00.000 2012-03-05 22:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2012-03-06 00:00:00.000 2012-03-06 00:00:01.000    SDO        AIA 171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    Length = 25 rows
    <BLANKLINE>
    <BLANKLINE>

To search for data from multiple instruments, wavelengths, times etc., use the pipe ``|`` operator which joins queries using a logical "OR" operator.
In this example we'll search for LYRA or RHESSI data in a given time range:

.. code-block:: python

    >>> Fido.search(a.Time('2012/3/4', '2012/3/4 02:00'),
    ...             a.Instrument.lyra | a.Instrument.rhessi)  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 3 Providers:
    <BLANKLINE>
    2 Results from the LYRAClient:
    Source: http://proba2.oma.be/lyra/data/bsd
    <BLANKLINE>
           Start Time               End Time        Instrument  Physobs   Source Provider Level
    ----------------------- ----------------------- ---------- ---------- ------ -------- -----
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999       LYRA irradiance PROBA2      ESA     2
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999       LYRA irradiance PROBA2      ESA     3
    <BLANKLINE>
    1 Results from the RHESSIClient:
    Source: https://hesperia.gsfc.nasa.gov/hessidata
    <BLANKLINE>
           Start Time               End Time        Instrument      Physobs       Source Provider
    ----------------------- ----------------------- ---------- ------------------ ------ --------
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999     RHESSI summary_lightcurve RHESSI     NASA
    <BLANKLINE>
    3 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        Source Instrument   Wavelength   Provider  Physobs  Extent Type   Size
                                                                           keV                                       Mibyte
    ----------------------- ----------------------- ------ ---------- -------------- -------- --------- ----------- --------
    2012-03-03 22:57:40.000 2012-03-04 00:33:20.000 RHESSI     RHESSI 3.0 .. 17000.0     LSSP intensity PARTIAL_SUN -0.00098
    2012-03-04 00:33:20.000 2012-03-04 01:45:40.000 RHESSI     RHESSI 3.0 .. 17000.0     LSSP intensity PARTIAL_SUN -0.00098
    2012-03-04 01:45:40.000 2012-03-04 02:09:00.000 RHESSI     RHESSI 3.0 .. 17000.0     LSSP intensity PARTIAL_SUN -0.00098
    <BLANKLINE>
    <BLANKLINE>


Working with Search Results
***************************

:meth:`Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>` can make multiple queries to multiple clients in one search.
This means that the results of a call to search can contain many sets of records, called responses, from many clients.
The results of a search are represented in a `~sunpy.net.fido_factory.UnifiedResponse` object, which provides access to all the response tables and allows some operations to be performed on all the results at once.
`~sunpy.net.fido_factory.UnifiedResponse` acts both like a two dimensional array, where the first dimension is the response index and the second index is the row index, and a dictionary where you can index the responses by the name of the client.

For example, the following code returns a response containing LYRA data from the `~sunpy.net.dataretriever.LYRAClient`, and EVE data from the `~sunpy.net.vso.VSOClient`:

.. code-block:: python

    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Level.two,
    ...                       a.Instrument.lyra | a.Instrument.eve)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    2 Results from the LYRAClient:
    Source: http://proba2.oma.be/lyra/data/bsd
    <BLANKLINE>
           Start Time               End Time        Instrument  Physobs   Source Provider Level
    ----------------------- ----------------------- ---------- ---------- ------ -------- -----
    2012-01-01 00:00:00.000 2012-01-01 23:59:59.999       LYRA irradiance PROBA2      ESA     2
    2012-01-02 00:00:00.000 2012-01-02 23:59:59.999       LYRA irradiance PROBA2      ESA     2
    <BLANKLINE>
    50 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        Source Instrument   Wavelength   Provider  Physobs   Extent Type   Size
                                                                         Angstrom                                     Mibyte
    ----------------------- ----------------------- ------ ---------- -------------- -------- ---------- ----------- --------
    2012-01-01 00:00:00.000 2012-01-01 01:00:00.000    SDO        EVE 93.0 .. 1033.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 00:00:00.000 2012-01-01 01:00:00.000    SDO        EVE 60.0 .. 1060.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 01:00:00.000 2012-01-01 02:00:00.000    SDO        EVE 93.0 .. 1033.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 01:00:00.000 2012-01-01 02:00:00.000    SDO        EVE 60.0 .. 1060.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 02:00:00.000 2012-01-01 03:00:00.000    SDO        EVE 93.0 .. 1033.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 02:00:00.000 2012-01-01 03:00:00.000    SDO        EVE 60.0 .. 1060.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 03:00:00.000 2012-01-01 04:00:00.000    SDO        EVE 93.0 .. 1033.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 03:00:00.000 2012-01-01 04:00:00.000    SDO        EVE 60.0 .. 1060.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 04:00:00.000 2012-01-01 05:00:00.000    SDO        EVE 93.0 .. 1033.0     LASP irradiance    FULLDISK -0.00098
                        ...                     ...    ...        ...            ...      ...        ...         ...      ...
    2012-01-01 20:00:00.000 2012-01-01 21:00:00.000    SDO        EVE 93.0 .. 1033.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 20:00:00.000 2012-01-01 21:00:00.000    SDO        EVE 60.0 .. 1060.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 21:00:00.000 2012-01-01 22:00:00.000    SDO        EVE 93.0 .. 1033.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 21:00:00.000 2012-01-01 22:00:00.000    SDO        EVE 60.0 .. 1060.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 22:00:00.000 2012-01-01 23:00:00.000    SDO        EVE 93.0 .. 1033.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 22:00:00.000 2012-01-01 23:00:00.000    SDO        EVE 60.0 .. 1060.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 23:00:00.000 2012-01-02 00:00:00.000    SDO        EVE 93.0 .. 1033.0     LASP irradiance    FULLDISK -0.00098
    2012-01-01 23:00:00.000 2012-01-02 00:00:00.000    SDO        EVE 60.0 .. 1060.0     LASP irradiance    FULLDISK -0.00098
    2012-01-02 00:00:00.000 2012-01-02 01:00:00.000    SDO        EVE 93.0 .. 1033.0     LASP irradiance    FULLDISK -0.00098
    2012-01-02 00:00:00.000 2012-01-02 01:00:00.000    SDO        EVE 60.0 .. 1060.0     LASP irradiance    FULLDISK -0.00098
    Length = 50 rows
    <BLANKLINE>
    <BLANKLINE>

If you then wanted to inspect just the LYRA data for the whole time range specified in the search, you would index this response to see just the results returned by the `~sunpy.net.dataretriever.LYRAClient`:

.. code-block:: python

    >>> results[0, :]  # doctest: +REMOTE_DATA
    <sunpy.net.dataretriever.client.QueryResponse object at ...>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-01-01 00:00:00.000 2012-01-01 23:59:59.999       LYRA ...      ESA     2
    2012-01-02 00:00:00.000 2012-01-02 23:59:59.999       LYRA ...      ESA     2

Or, equivalently:

.. code-block:: python

    >>> results["lyra"]  # doctest: +REMOTE_DATA
    <sunpy.net.dataretriever.client.QueryResponse object at ...>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-01-01 00:00:00.000 2012-01-01 23:59:59.999       LYRA ...      ESA     2
    2012-01-02 00:00:00.000 2012-01-02 23:59:59.999       LYRA ...      ESA     2

Normal slicing operations work as with any other Python sequence, e.g. ``results[1, ::10]`` to access every tenth file in the result returned by the second client.

Note that the first (response) index is still necessary even if results are only found for a single client.
So in this case the first result would be ``results[0, 0]`` rather than ``results[0]`` (the latter would return all results from the first - and only - client and is therefore the same as ``results``).

As we have seen above the `~sunpy.net.fido_factory.UnifiedResponse` object contains many response tables which make up the search results.
Each of the responses are `~sunpy.net.base_client.QueryResponseTable` objects, which are `astropy.table` objects meaning that you can interact with them and filter them like any other tabular data.
This can be used to interact with results which are metadata only, i.e. searches from the HEK, or it can be used to reduce the number of files downloaded by `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`.

For example if we did a query for some AIA and HMI data:

.. code-block:: python

    >>> results = Fido.search(a.Time("2020/01/01", "2020/01/01 00:01"), a.Instrument.aia | a.Instrument.hmi)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    41 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 2.779 Gbyte
    <BLANKLINE>
           Start Time               End Time        Source Instrument    Wavelength    Provider  Physobs  Wavetype Extent Width Extent Length Extent Type   Size
                                                                          Angstrom                                                                         Mibyte
    ----------------------- ----------------------- ------ ---------- ---------------- -------- --------- -------- ------------ ------------- ----------- --------
    2020-01-01 00:00:00.000 2020-01-01 00:00:01.000    SDO        AIA   335.0 .. 335.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:04.000 2020-01-01 00:00:05.000    SDO        AIA   193.0 .. 193.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:05.000 2020-01-01 00:00:06.000    SDO        AIA   304.0 .. 304.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:05.000 2020-01-01 00:00:06.000    SDO        AIA 4500.0 .. 4500.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:06.000 2020-01-01 00:00:07.000    SDO        AIA   131.0 .. 131.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:09.000 2020-01-01 00:00:10.000    SDO        AIA   171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:09.000 2020-01-01 00:00:10.000    SDO        AIA   211.0 .. 211.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:11.000 2020-01-01 00:00:12.000    SDO        AIA     94.0 .. 94.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:12.000 2020-01-01 00:00:13.000    SDO        AIA   335.0 .. 335.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
                        ...                     ...    ...        ...              ...      ...       ...      ...          ...           ...         ...      ...
    2020-01-01 00:00:47.000 2020-01-01 00:00:48.000    SDO        AIA     94.0 .. 94.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:48.000 2020-01-01 00:00:49.000    SDO        AIA   335.0 .. 335.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:52.000 2020-01-01 00:00:53.000    SDO        AIA 1700.0 .. 1700.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:52.000 2020-01-01 00:00:53.000    SDO        AIA   193.0 .. 193.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:53.000 2020-01-01 00:00:54.000    SDO        AIA   304.0 .. 304.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:54.000 2020-01-01 00:00:55.000    SDO        AIA   131.0 .. 131.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:57.000 2020-01-01 00:00:58.000    SDO        AIA   171.0 .. 171.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:57.000 2020-01-01 00:00:58.000    SDO        AIA   211.0 .. 211.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:00:59.000 2020-01-01 00:01:00.000    SDO        AIA     94.0 .. 94.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    2020-01-01 00:01:00.000 2020-01-01 00:01:01.000    SDO        AIA   335.0 .. 335.0     JSOC intensity   NARROW         4096          4096    FULLDISK 64.64844
    Length = 41 rows
    <BLANKLINE>
    3 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        Source Instrument    Wavelength    Provider      Physobs       Wavetype Extent Width Extent Length Extent Type   Size
                                                                          Angstrom                                                                                  Mibyte
    ----------------------- ----------------------- ------ ---------- ---------------- -------- ------------------ -------- ------------ ------------- ----------- --------
    2020-01-01 00:00:22.000 2020-01-01 00:00:23.000    SDO        HMI 6173.0 .. 6174.0     JSOC          intensity   NARROW         4096          4096    FULLDISK -0.00098
    2020-01-01 00:00:22.000 2020-01-01 00:00:23.000    SDO        HMI 6173.0 .. 6174.0     JSOC LOS_magnetic_field   NARROW         4096          4096    FULLDISK -0.00098
    2020-01-01 00:00:22.000 2020-01-01 00:00:23.000    SDO        HMI 6173.0 .. 6174.0     JSOC       LOS_velocity   NARROW         4096          4096    FULLDISK -0.00098
    <BLANKLINE>
    <BLANKLINE>

The VSO client returns a lot of information about the records, so the first thing we can do is show only the columns we are interested in.
We can inspect all the available column names in all the responses with the `~.UnifiedResponse.all_colnames` property:

.. code-block:: python

    >>> results.all_colnames  # doctest: +REMOTE_DATA
    ['End Time', 'Extent Length', 'Extent Type', 'Extent Width', 'Instrument', 'Physobs', 'Provider', 'Size', 'Source', 'Start Time', 'Wavelength', 'Wavetype', 'fileid']

And then we can pick which ones to see with the :meth:`~.UnifiedResponse.show` method:

.. code-block:: python

    >>> results.show("Start Time", "Instrument", "Physobs", "Wavelength")  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    41 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time       Instrument  Physobs     Wavelength
                                                     Angstrom
    ----------------------- ---------- --------- ----------------
    2020-01-01 00:00:00.000        AIA intensity   335.0 .. 335.0
    2020-01-01 00:00:04.000        AIA intensity   193.0 .. 193.0
    2020-01-01 00:00:05.000        AIA intensity   304.0 .. 304.0
    2020-01-01 00:00:05.000        AIA intensity 4500.0 .. 4500.0
    2020-01-01 00:00:06.000        AIA intensity   131.0 .. 131.0
    2020-01-01 00:00:09.000        AIA intensity   171.0 .. 171.0
    2020-01-01 00:00:09.000        AIA intensity   211.0 .. 211.0
    2020-01-01 00:00:11.000        AIA intensity     94.0 .. 94.0
    2020-01-01 00:00:12.000        AIA intensity   335.0 .. 335.0
                        ...        ...       ...              ...
    2020-01-01 00:00:47.000        AIA intensity     94.0 .. 94.0
    2020-01-01 00:00:48.000        AIA intensity   335.0 .. 335.0
    2020-01-01 00:00:52.000        AIA intensity 1700.0 .. 1700.0
    2020-01-01 00:00:52.000        AIA intensity   193.0 .. 193.0
    2020-01-01 00:00:53.000        AIA intensity   304.0 .. 304.0
    2020-01-01 00:00:54.000        AIA intensity   131.0 .. 131.0
    2020-01-01 00:00:57.000        AIA intensity   171.0 .. 171.0
    2020-01-01 00:00:57.000        AIA intensity   211.0 .. 211.0
    2020-01-01 00:00:59.000        AIA intensity     94.0 .. 94.0
    2020-01-01 00:01:00.000        AIA intensity   335.0 .. 335.0
    Length = 41 rows
    <BLANKLINE>
    3 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time       Instrument      Physobs          Wavelength
                                                              Angstrom
    ----------------------- ---------- ------------------ ----------------
    2020-01-01 00:00:22.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:00:22.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:00:22.000        HMI       LOS_velocity 6173.0 .. 6174.0
    <BLANKLINE>
    <BLANKLINE>

To give an example of filtering post-search, let's only return the rows in the table which are line-of-sight magnetograms from HMI or the 94Å passband from AIA.
You can also always do this filtering with the `a.Physobs <sunpy.net.attrs.Physobs>` and `a.Wavelength <sunpy.net.attrs.Wavelength>` attrs in the search command.

First we split the results in to a table for AIA and a table for HMI:

.. code-block:: python

   >>> aia, hmi = results  # doctest: +REMOTE_DATA

We can use boolean indexing to match the value of the ``"Physobs"`` column:

.. code-block:: python

    >>> hmi_los = hmi[hmi["Physobs"] == "LOS_magnetic_field"]  # doctest: +REMOTE_DATA
    >>> hmi_los.show("Start Time", "Instrument", "Wavelength", "Physobs")  # doctest: +REMOTE_DATA
    <sunpy.net.vso.table_response.VSOQueryResponseTable object at ...>
           Start Time       Instrument    Wavelength         Physobs
                                           Angstrom
    ----------------------- ---------- ---------------- ------------------
    2020-01-01 00:00:22.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field

To match the ``"Wavelength"`` column we need to account for the fact that VSO results return a wavelength range of ``[min, max]`` so we match the min:

.. code-block:: python

    >>> aia_94 = aia[aia["Wavelength"][:, 0] == 94 * u.AA]  # doctest: +REMOTE_DATA
    >>> aia_94.show("Start Time", "Instrument", "Wavelength", "Physobs")  # doctest: +REMOTE_DATA
    <sunpy.net.vso.table_response.VSOQueryResponseTable object at ...>
           Start Time       Instrument  Wavelength   Physobs
                                         Angstrom
    ----------------------- ---------- ------------ ---------
    2020-01-01 00:00:11.000        AIA 94.0 .. 94.0 intensity
    2020-01-01 00:00:23.000        AIA 94.0 .. 94.0 intensity
    2020-01-01 00:00:35.000        AIA 94.0 .. 94.0 intensity
    2020-01-01 00:00:47.000        AIA 94.0 .. 94.0 intensity
    2020-01-01 00:00:59.000        AIA 94.0 .. 94.0 intensity

.. warning::

   While you can reduce the number of columns and rows in the results, the ``fetch()`` method that downloads data may need certain columns to be present to successfully download the files.
   It is therefore highly recommended that if you are planning on downloading data you do not slice out columns, but instead use ``.show()`` to only display the ones you are interested in.

.. _sunpy-tutorial-acquiring-data-downloading-data:

Downloading data
****************

Once you have located your files via a `Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>`, you can download them via `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`.
Here we'll just download the first file in the result:

.. code-block:: python

    >>> downloaded_files = Fido.fetch(results[0, 0]) # doctest: +REMOTE_DATA
    >>> downloaded_files # doctest: +REMOTE_DATA
    <parfive.results.Results object at ...>
    ['.../aia.lev1.335A_2020_01_01T00_00_00.64Z.image_lev1.fits']

This downloads the files to the location set in the sunpy config file.
It also returns a `parfive.Results` object ``downloaded_files``, which contains local file paths to all the downloaded data.

You can also explicitly specify the path to which you want the data downloaded:

.. code-block:: python

  >>> downloaded_files = Fido.fetch(results, path='/ThisIs/MyPath/to/Data/{file}')  # doctest: +SKIP

This downloads the query results into the directory ``/ThisIs/MyPath/to/Data``, naming each downloaded file with the filename ``{file}`` obtained from the client.
You can also use other properties of the returned query to define the path where the data is saved.
For example, to save the data to a subdirectory named after the instrument, use:

.. code-block:: python

    >>> downloaded_files = Fido.fetch(results, path='./{instrument}/{file}')  # doctest: +SKIP

You can see the list of options that can be specified in path for all the files to be downloaded with ``results.path_format_keys``:

.. code-block:: python

    >>> sorted(results.path_format_keys()) # doctest: +REMOTE_DATA
    ['end_time', 'extent_length', 'extent_type', 'extent_width', 'fileid', 'instrument', 'physobs', 'provider', 'size', 'source', 'start_time', 'wavelength', 'wavetype']

Retrying Downloads
^^^^^^^^^^^^^^^^^^

If any files failed to download, the progress bar will show an incomplete number of files (i.e. 100/150) and the `parfive.Results` object will contain a list of the URLs that failed to transfer and the error associated with them.
This can be accessed with the ``.errors`` attribute or by printing the `~parfive.Results` object:

.. code-block:: python

    >>> print(downloaded_files.errors)  # doctest: +SKIP

The transfer can be retried by passing the `parfive.Results` object back to `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`:

.. code-block:: python

    >>> downloaded_files = Fido.fetch(downloaded_files)  # doctest: +SKIP

doing this will append any newly downloaded file names to the list and replace the ``.errors`` list with any errors that occurred during the second attempt.
