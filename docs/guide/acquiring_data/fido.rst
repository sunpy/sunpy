.. _fido_guide:

***************************************
Finding and Downloading Data using Fido
***************************************

This guide outlines how to search for and download data using SunPy's
Federated Internet Data Obtainer...or more usually (and sanely) referred to as Fido.
`~sunpy.net.Fido` is a unified interface for searching
and fetching solar physics data irrespective of the underlying
client or webservice through which the data is obtained, e.g. VSO_,
JSOC_, etc. It therefore supplies a single, easy and consistent way to
obtain most forms of solar physics data.

Import
******

The `~sunpy.net.Fido` object is in `sunpy.net`.
It can be imported as follows::

    >>> from sunpy.net import Fido, attrs as a

Search Attributes
*****************

To search for data with `~sunpy.net.Fido`, you need to specify attributes to search against.
The range of attributes are found in the `attrs <sunpy.net.attrs>` submodule.
Examples of these attributes are:

- `a.Time <sunpy.net.attrs.Time>`
- `a.Instrument <sunpy.net.attrs.Instrument>`
- `a.Wavelength <sunpy.net.attrs.Wavelength>`

whereas some of these attributes are client specific, and are found under `attrs.vso <sunpy.net.vso.attrs>` and `attrs.jsoc <sunpy.net.jsoc.attrs>`.

In to each attribute you have to provide a value to use::

    >>> a.Time('2012/3/4', '2012/3/6'), a.Instrument.lyra
    (<sunpy.net.attrs.Time(2012-03-04 00:00:00.000, 2012-03-06 00:00:00.000)>, <sunpy.net.attrs.Instrument(LYRA: Lyman Alpha Radiometer is the solar UV radiometer on board
    Proba-2.) object at ...>)

For attributes that have no fixed selection of values (``Time`` for example) you will have to provide the range you require.
However, for attributes that have a fixed range of **known** values, it is possible to list all these values.
**Please note that each list is not exhaustive.**

Using ``Instrument`` as the first example, if you print the object::

    >>> print(a.Instrument)
    sunpy.net.attrs.Instrument
    <BLANKLINE>
    Specifies the Instrument name for the search.
    <BLANKLINE>
           Attribute Name          Client          Full Name                                           Description
    --------------------------- ----------- ------------------------ --------------------------------------------------------------------------------
    aia                         VSO         AIA                      Atmospheric Imaging Assembly
    bbi                         VSO         BBI                      None
    bcs                         VSO         BCS                      Bragg Crystal Spectrometer
    be_continuum                VSO         BE-Continuum             INAF-OACT Barra Equatoriale Continuum Instrument
    be_halpha                   VSO         BE-Halpha                INAF-OACT Barra Equatoriale Hα Instrument
    bic_hifi                    VSO         BIC-HIFI                 None
    bigbear                     VSO         Big Bear                 Big Bear Solar Observatory, California TON and GONG+ sites
    caii                        VSO         CAII                     Kanzelhöhe Ca II k Instrument
    cds                         VSO         CDS                      Coronal Diagnostic Spectrometer
    ...

You get a full list of known values, a description and what "Clients" support those values (if you want to use a specific data source).
This is supported for most attributes including the client specific ones.


For JSOC::

    >>> print(a.jsoc.Series)
    sunpy.net.jsoc.attrs.Series
    <BLANKLINE>
    The JSOC Series to Download.
    <BLANKLINE>
              Attribute Name           Client             Full Name                                                Description
    ---------------------------------- ------ ---------------------------------- --------------------------------------------------------------------------------
    aia_flatfield                      JSOC   aia.flatfield                      AIA flatfield
    aia_lev1                           JSOC   aia.lev1                           AIA Level 1
    aia_lev1_euv_12s                   JSOC   aia.lev1_euv_12s                   AIA Level 1, 12 second cadence
    aia_lev1_uv_24s                    JSOC   aia.lev1_uv_24s                    AIA Level 1, 24 second cadence
    aia_lev1_vis_1h                    JSOC   aia.lev1_vis_1h                    AIA Level 1, 3600 second cadence
    aia_master_pointing3h              JSOC   aia.master_pointing3h              Master Pointing Parameters
    aia_response                       JSOC   aia.response                       AIA instrument response table
    aia_temperature_summary_300s       JSOC   aia.temperature_summary_300s       Temperature Statistics from AIA Housekeeping - Thermal Packet
    hmi_b_135s                         JSOC   hmi.b_135s                         Full-disk Milne-Eddington inversion with the azimuth disambiguation informati...
    ...

Furthermore, you can use tab completion to auto-fill the attribute name, for example by typing ``a.jsoc.aia_f<TAB>``.

Searching for Data Using Fido
*****************************

For example::

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.lyra, a.Level.two) # doctest: +REMOTE_DATA

this returns an `~sunpy.net.fido_factory.UnifiedResponse` object containing
information on the available online files which fit the criteria specified by
the attrs objects in the above call. It does not download the files. For
instructions on how to download data using Fido, see :ref:`downloading_data`.

To see a summary of results of our query, simple type the name of the
variable set to the Fido search, in this case, result::

    >>> result  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the LYRAClient:
         Start Time           End Time      Instrument ... Source Provider Level
    ------------------- ------------------- ---------- ... ------ -------- -----
    2012-03-04 00:00:00 2012-03-04 23:59:59       LYRA ... PROBA2      ESA     2
    2012-03-05 00:00:00 2012-03-05 23:59:59       LYRA ... PROBA2      ESA     2
    2012-03-06 00:00:00 2012-03-06 23:59:59       LYRA ... PROBA2      ESA     2
    <BLANKLINE>
    <BLANKLINE>

Queries can be made more flexible or specific by adding more attrs objects to
the `~sunpy.net.Fido` search. Specific
passbands can be searched for by supplying an `~astropy.units.Quantity` to the
`a.Wavelength <sunpy.net.attrs.Wavelength>` attribute::

    >>> import astropy.units as u
    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.norh,
    ...             a.Wavelength(17*u.GHz))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the NoRHClient:
         Start Time           End Time      Instrument Source Provider Wavelength
    ------------------- ------------------- ---------- ------ -------- ----------
    2012-03-04 00:00:00 2012-03-04 23:59:59       NORH   NAOJ      NRO   17.0 GHz
    2012-03-05 00:00:00 2012-03-05 23:59:59       NORH   NAOJ      NRO   17.0 GHz
    2012-03-06 00:00:00 2012-03-06 23:59:59       NORH   NAOJ      NRO   17.0 GHz
    <BLANKLINE>
    <BLANKLINE>

Data of a given cadence can also be specified using the Sample attribute. To
search for data at a given cadence use the
`a.Sample <sunpy.net.attrs.Sample>` attribute.
`a.Sample <sunpy.net.attrs.Sample>` is only supported by the
`sunpy.net.vso.VSOClient` hence it has the ``a.vso`` prefix. Attributes
like this which are client specific will result in
`~sunpy.net.Fido` only searching that
client for results, in this case VSO.::

    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.aia,
    ...             a.Wavelength(171*u.angstrom), a.Sample(10*u.minute))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    289 Results from the VSOClient:
        Start Time [1]       End Time [1]    Source ...   Type   Wavelength [2]
                                                    ...             Angstrom
     ------------------- ------------------- ------ ... -------- --------------
     2012-03-04 00:00:00 2012-03-04 00:00:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-04 00:10:00 2012-03-04 00:10:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-04 00:20:00 2012-03-04 00:20:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-04 00:30:00 2012-03-04 00:30:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-04 00:40:00 2012-03-04 00:40:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-04 00:50:00 2012-03-04 00:50:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-04 01:00:00 2012-03-04 01:00:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-04 01:10:00 2012-03-04 01:10:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-04 01:20:00 2012-03-04 01:20:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-04 01:30:00 2012-03-04 01:30:01    SDO ... FULLDISK 171.0 .. 171.0
                     ...                 ...    ... ...      ...            ...
     2012-03-05 22:30:00 2012-03-05 22:30:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-05 22:40:00 2012-03-05 22:40:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-05 22:50:00 2012-03-05 22:50:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-05 23:00:00 2012-03-05 23:00:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-05 23:10:00 2012-03-05 23:10:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-05 23:20:00 2012-03-05 23:20:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-05 23:30:00 2012-03-05 23:30:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-05 23:40:00 2012-03-05 23:40:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-05 23:50:00 2012-03-05 23:50:01    SDO ... FULLDISK 171.0 .. 171.0
     2012-03-06 00:00:00 2012-03-06 00:00:01    SDO ... FULLDISK 171.0 .. 171.0
    Length = 289 rows
    <BLANKLINE>
    <BLANKLINE>

To search for data from multiple instruments, wavelengths, times etc., use the
pipe ``|`` operator. This joins queries together just as the logical ``OR``
operator would::

    >>> Fido.search(a.Time('2012/3/4', '2012/3/4 02:00'),
    ...             a.Instrument.lyra | a.Instrument.rhessi)  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 3 Providers:
    <BLANKLINE>
    2 Results from the LYRAClient:
         Start Time           End Time      Instrument ... Source Provider Level
    ------------------- ------------------- ---------- ... ------ -------- -----
    2012-03-04 00:00:00 2012-03-04 23:59:59       LYRA ... PROBA2      ESA     2
    2012-03-04 00:00:00 2012-03-04 23:59:59       LYRA ... PROBA2      ESA     3
    <BLANKLINE>
    1 Results from the RHESSIClient:
         Start Time           End Time      Instrument ... Source Provider
    ------------------- ------------------- ---------- ... ------ --------
    2012-03-04 00:00:00 2012-03-04 23:59:59     RHESSI ... RHESSI     NASA
    <BLANKLINE>
    3 Results from the VSOClient:
       Start Time [1]       End Time [1]    Source ...     Type    Wavelength [2]
                                                   ...                  keV
    ------------------- ------------------- ------ ... ----------- --------------
    2012-03-03 22:57:40 2012-03-04 00:33:20 RHESSI ... PARTIAL_SUN 3.0 .. 17000.0
    2012-03-04 00:33:20 2012-03-04 01:45:40 RHESSI ... PARTIAL_SUN 3.0 .. 17000.0
    2012-03-04 01:45:40 2012-03-04 02:09:00 RHESSI ... PARTIAL_SUN 3.0 .. 17000.0
    <BLANKLINE>
    <BLANKLINE>

Indexing search results
***********************

The `~sunpy.net.fido_factory.UnifiedResponse` that Fido returns can be
indexed to access a subset of the search results. When doing this, the
results should be treated as a two-dimensional array in which the first
dimension corresponds to the clients which have returned results and the
second to the records returned.

For example, the following code returns a response containing LYRA data from
the `~sunpy.net.dataretriever.LYRAClient`, and EVE data from the
`~sunpy.net.vso.VSOClient`::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Level.two,
    ...                       a.Instrument.lyra | a.Instrument.eve)  # doctest: +REMOTE_DATA

If you then wanted to inspect just the LYRA data for the whole time range
specified in the search, you would index this response to see just the
results returned by the `~sunpy.net.dataretriever.LYRAClient`::

    >>> results[0, :]  # doctest: +REMOTE_DATA
    <sunpy.net.dataretriever.client.QueryResponse object at ...>
         Start Time           End Time      Instrument ... Source Provider Level
    ------------------- ------------------- ---------- ... ------ -------- -----
    2012-01-01 00:00:00 2012-01-01 23:59:59       LYRA ... PROBA2      ESA     2
    2012-01-02 00:00:00 2012-01-02 23:59:59       LYRA ... PROBA2      ESA     2

Or, equivalently::

    >>> results[0]  # doctest: +REMOTE_DATA
    <sunpy.net.dataretriever.client.QueryResponse object at ...>
         Start Time           End Time      Instrument ... Source Provider Level
    ------------------- ------------------- ---------- ... ------ -------- -----
    2012-01-01 00:00:00 2012-01-01 23:59:59       LYRA ... PROBA2      ESA     2
    2012-01-02 00:00:00 2012-01-02 23:59:59       LYRA ... PROBA2      ESA     2

Normal slicing operations work as with any other Python sequence, e.g.
``results[1,::10]`` to access every tenth file in the result returned by
the second client.

Note that the first (client) index is still necessary even if results
are only found for a single client. So in this case the first result
would be ``results[0,0]`` rather than ``results[0]`` (the latter would return
all results from the first - and only - client and is therefore the
same as ``results``).

.. _downloading_data:

Downloading data
****************
Once you have located your files via a
`Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>`, you can
download them via `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`::

    >>> downloaded_files = Fido.fetch(results)  # doctest: +SKIP

This downloads the files to the location set in you sunpy config file. It also
returns a `parfive.Results` object ``downloaded_files``, of absolute file paths
of where the files have been downloaded to.

You can also specify the path to which you want the data downloaded::

  >>> downloaded_files = Fido.fetch(results, path='/ThisIs/MyPath/to/Data/{file}')  # doctest: +SKIP

This downloads the query results into the directory
``/ThisIs/MyPath/to/Data``, naming each downloaded file with the
filename ``{file}`` obtained from the client.
You can also use other properties of the returned query
to define the path where the data is saved.  For example, to save the
data to a subdirectory named after the instrument, use::

    >>> downloaded_files = Fido.fetch(results, path='./{instrument}/{file}')  # doctest: +SKIP

You can see the list of options that can be specified in path for all the files
to be downloaded with ``results.response_block_properties``.

Retrying Downloads
==================

If any files failed to download, the progress bar will show an incomplete number
of files (i.e. 100/150) and the `parfive.Results` object will contain a list of
the URLs that failed to transfer and the error associated with them. This can be
accessed with the ``.errors`` attribute or by printing the `~parfive.Results`
object::

    >>> print(downloaded_files.errors)  # doctest: +SKIP

The transfer can be retried by passing the `parfive.Results` object back to
`Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`::

    >>> downloaded_files = Fido.fetch(downloaded_files)  # doctest: +SKIP

doing this will append any newly downloaded file names to the list and replace
the ``.errors`` list with any errors that occurred during the second attempt.


.. _VSO: https://sdac.virtualsolar.org/cgi/search
.. _JSOC: http://jsoc.stanford.edu/


Fido Clients
************

`~sunpy.net.Fido` provides access to many sources of data via "clients", these clients can be defined inside sunpy or in other packages.
If you want to see the current list of clients you can do::

    >>> print(Fido)
    sunpy.net.Fido
    <BLANKLINE>
    Fido is a unified data search and retrieval tool.
    <BLANKLINE>
    It provides simultaneous access to a variety of online data sources, some
    cover multiple instruments and data products like the Virtual Solar
    Observatory and some are specific to a single source.
    <BLANKLINE>
    For details of using `~sunpy.net.Fido` see :ref:`fido_guide`.
    <BLANKLINE>
    <BLANKLINE>
          Client                                                    Description
    ----------------- -------------------------------------------------------------------------------------------------------
    EVEClient         Provides access to Level 0C Extreme ultraviolet Variability Experiment (EVE) data.
    GBMClient         Provides access to data from the Gamma-Ray Burst Monitor (GBM) instrument on board the Fermi satellite.
    XRSClient         Provides access to the GOES XRS fits files archive.
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
