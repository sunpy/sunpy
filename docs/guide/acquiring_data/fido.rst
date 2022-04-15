.. _fido_guide:

***************************************
Finding and Downloading Data using Fido
***************************************

This guide outlines how to search for and download data using `~sunpy.net.Fido` sunpy's interface for search and download.
`~sunpy.net.Fido` is a unified interface for searching and fetching solar physics data irrespective of the underlying client or webservice through which the data is obtained, e.g. VSO_,JSOC_, etc.
It therefore supplies a single, easy and consistent way to obtain most forms of solar physics data.

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

whereas some of these attributes are client specific, and are found under client specific submodules, e.g. `attrs.vso <sunpy.net.vso.attrs>` and `attrs.jsoc <sunpy.net.jsoc.attrs>`.

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
    bcs                         VSO         BCS                      Bragg Crystal Spectrometer
    be_continuum                VSO         BE-Continuum             INAF-OACT Barra Equatoriale Continuum Instrument
    be_halpha                   VSO         BE-Halpha                INAF-OACT Barra Equatoriale Hα Instrument
    bigbear                     VSO         Big Bear                 Big Bear Solar Observatory, California TON and GONG+ sites
    caii                        VSO         CAII                     Kanzelhöhe Ca II k Instrument
    cds                         VSO         CDS                      Coronal Diagnostic Spectrometer
    celias                      VSO         CELIAS                   Charge, Element, and Isotope Analysis System
    ...

You get a full list of known values, a description and what "Clients" support those values (if you want to use a specific data source).
This is supported for most attributes including the client specific ones.


For example you can print a list of Series provided by JSOC::

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

this returns an `~sunpy.net.fido_factory.UnifiedResponse` object containing information on the results which fit the criteria specified by the attrs objects in the above call.
It does not download the files.
For instructions on how to download data using Fido, see :ref:`downloading_data`.

To see a summary of the results of our query, simply type the name of the variable set to the Fido search, in this case, result::

    >>> result  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
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
Specific passbands can be searched for by supplying an `~astropy.units.Quantity` to the `a.Wavelength <sunpy.net.attrs.Wavelength>` attribute::

    >>> import astropy.units as u
    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.norh,
    ...             a.Wavelength(17*u.GHz))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    3 Results from the NoRHClient:
    Source: https://solar.nro.nao.ac.jp/norh/doc/manuale/node1.html
    <BLANKLINE>
           Start Time               End Time        ... Provider Wavelength
                                                    ...             GHz
    ----------------------- ----------------------- ... -------- ----------
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999 ...      NRO       17.0
    2012-03-05 00:00:00.000 2012-03-05 23:59:59.999 ...      NRO       17.0
    2012-03-06 00:00:00.000 2012-03-06 23:59:59.999 ...      NRO       17.0
    <BLANKLINE>
    <BLANKLINE>

Data of a given cadence can also be specified using the Sample attribute.
To search for data at a given cadence use the `a.Sample <sunpy.net.attrs.Sample>` attribute.

    >>> Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument.aia,
    ...             a.Wavelength(171*u.angstrom), a.Sample(10*u.minute))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    289 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 19.591 Gbyte
    <BLANKLINE>
           Start Time       ...
                            ...
    ----------------------- ...
    2012-03-04 00:00:00.000 ...
    2012-03-04 00:10:00.000 ...
    2012-03-04 00:20:00.000 ...
    2012-03-04 00:30:00.000 ...
    2012-03-04 00:40:00.000 ...
    2012-03-04 00:50:00.000 ...
    2012-03-04 01:00:00.000 ...
    2012-03-04 01:10:00.000 ...
    2012-03-04 01:20:00.000 ...
    2012-03-04 01:30:00.000 ...
                        ... ...
    2012-03-05 22:30:00.000 ...
    2012-03-05 22:40:00.000 ...
    2012-03-05 22:50:00.000 ...
    2012-03-05 23:00:00.000 ...
    2012-03-05 23:10:00.000 ...
    2012-03-05 23:20:00.000 ...
    2012-03-05 23:30:00.000 ...
    2012-03-05 23:40:00.000 ...
    2012-03-05 23:50:00.000 ...
    2012-03-06 00:00:00.000 ...
    Length = 289 rows
    <BLANKLINE>
    <BLANKLINE>

To search for data from multiple instruments, wavelengths, times etc., use the pipe ``|`` operator.
This joins queries together just as the logical ``OR`` operator would::

    >>> Fido.search(a.Time('2012/3/4', '2012/3/4 02:00'),
    ...             a.Instrument.lyra | a.Instrument.rhessi)  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 3 Providers:
    <BLANKLINE>
    2 Results from the LYRAClient:
    Source: http://proba2.oma.be/lyra/data/bsd
    <BLANKLINE>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999       LYRA ...      ESA     2
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999       LYRA ...      ESA     3
    <BLANKLINE>
    1 Results from the RHESSIClient:
    Source: https://hesperia.gsfc.nasa.gov/hessidata
    <BLANKLINE>
           Start Time               End Time        Instrument ... Source Provider
    ----------------------- ----------------------- ---------- ... ------ --------
    2012-03-04 00:00:00.000 2012-03-04 23:59:59.999     RHESSI ... RHESSI     NASA
    <BLANKLINE>
    3 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        ...   Size        Info
                                                    ...  Mibyte
    ----------------------- ----------------------- ... -------- --------------
    2012-03-03 22:57:40.000 2012-03-04 00:33:20.000 ... -0.00098 RHESSI level-0
    2012-03-04 00:33:20.000 2012-03-04 01:45:40.000 ... -0.00098 RHESSI level-0
    2012-03-04 01:45:40.000 2012-03-04 02:09:00.000 ... -0.00098 RHESSI level-0
    <BLANKLINE>
    <BLANKLINE>


Working with Search Results
***************************

:meth:`Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>` can make multiple queries to multiple clients in one search.
This means that the results of a call to search can contain many sets of records, called responses, from many clients.
The results of a search are represented in a `~sunpy.net.fido_factory.UnifiedResponse` object, which provides access to all the response tables and allows some operations to be performed on all the results at once.
`~sunpy.net.fido_factory.UnifiedResponse` acts both like a two dimensional array, where the first dimension is the response index and the second index is the row index, and a dictionary where you can index the responses by the name of the client.

For example, the following code returns a response containing LYRA data from the `~sunpy.net.dataretriever.LYRAClient`, and EVE data from the `~sunpy.net.vso.VSOClient`::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Level.two,
    ...                       a.Instrument.lyra | a.Instrument.eve)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    2 Results from the LYRAClient:
    Source: http://proba2.oma.be/lyra/data/bsd
    <BLANKLINE>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-01-01 00:00:00.000 2012-01-01 23:59:59.999       LYRA ...      ESA     2
    2012-01-02 00:00:00.000 2012-01-02 23:59:59.999       LYRA ...      ESA     2
    <BLANKLINE>
    50 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        ...   Size         Info
                                                    ...  Mibyte
    ----------------------- ----------------------- ... -------- ----------------
    2012-01-01 00:00:00.000 2012-01-01 01:00:00.000 ... -0.00098 L2Lines (merged)
    2012-01-01 00:00:00.000 2012-01-01 01:00:00.000 ... -0.00098 L2Spectra (MEGS)
    2012-01-01 01:00:00.000 2012-01-01 02:00:00.000 ... -0.00098 L2Lines (merged)
                        ...                     ... ...      ...              ...
    2012-01-01 23:00:00.000 2012-01-02 00:00:00.000 ... -0.00098 L2Spectra (MEGS)
    2012-01-02 00:00:00.000 2012-01-02 01:00:00.000 ... -0.00098 L2Lines (merged)
    2012-01-02 00:00:00.000 2012-01-02 01:00:00.000 ... -0.00098 L2Spectra (MEGS)
    Length = 50 rows
    <BLANKLINE>
    <BLANKLINE>


If you then wanted to inspect just the LYRA data for the whole time range specified in the search, you would index this response to see just the results returned by the `~sunpy.net.dataretriever.LYRAClient`::

    >>> results[0, :]  # doctest: +REMOTE_DATA
    <sunpy.net.dataretriever.client.QueryResponse object at ...>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-01-01 00:00:00.000 2012-01-01 23:59:59.999       LYRA ...      ESA     2
    2012-01-02 00:00:00.000 2012-01-02 23:59:59.999       LYRA ...      ESA     2

Or, equivalently::

    >>> results["lyra"]  # doctest: +REMOTE_DATA
    <sunpy.net.dataretriever.client.QueryResponse object at ...>
           Start Time               End Time        Instrument ... Provider Level
    ----------------------- ----------------------- ---------- ... -------- -----
    2012-01-01 00:00:00.000 2012-01-01 23:59:59.999       LYRA ...      ESA     2
    2012-01-02 00:00:00.000 2012-01-02 23:59:59.999       LYRA ...      ESA     2

Normal slicing operations work as with any other Python sequence, e.g. ``results[1,::10]`` to access every tenth file in the result returned by the second client.

Note that the first (response) index is still necessary even if results are only found for a single client.
So in this case the first result would be ``results[0, 0]`` rather than ``results[0]`` (the latter would return all results from the first - and only - client and is therefore the same as ``results``).

Working with Response Tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As we have seen above the `~sunpy.net.fido_factory.UnifiedResponse` object contains many response tables which make up the search results.
Each of the responses are `~sunpy.net.base_client.QueryResponseTable` objects, which are `astropy.table` objects meaning that you can interact with them and filter them like any other tabular data.
This can be used to interact with results which are metadata only, i.e. searches from the HEK, or it can be used to reduce the number of files downloaded by `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`.

For example if we did a query for some AIA and HMI data::

    >>> results = Fido.search(a.Time("2020/01/01", "2020/01/01 00:05"), a.Instrument.aia | a.Instrument.hmi)  # doctest: +REMOTE_DATA
    >>> results  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    201 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 13.626 Gbyte
    <BLANKLINE>
           Start Time       ...
                            ...
    ----------------------- ...
    2020-01-01 00:00:00.000 ...
    2020-01-01 00:00:04.000 ...
    2020-01-01 00:00:05.000 ...
    2020-01-01 00:00:05.000 ...
    2020-01-01 00:00:06.000 ...
    2020-01-01 00:00:09.000 ...
    2020-01-01 00:00:09.000 ...
    2020-01-01 00:00:11.000 ...
    2020-01-01 00:00:12.000 ...
    2020-01-01 00:00:14.000 ...
                        ... ...
    2020-01-01 00:04:47.000 ...
    2020-01-01 00:04:48.000 ...
    2020-01-01 00:04:52.000 ...
    2020-01-01 00:04:52.000 ...
    2020-01-01 00:04:53.000 ...
    2020-01-01 00:04:54.000 ...
    2020-01-01 00:04:57.000 ...
    2020-01-01 00:04:57.000 ...
    2020-01-01 00:04:59.000 ...
    2020-01-01 00:05:00.000 ...
    Length = 201 rows
    <BLANKLINE>
    21 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        ...            Info
                                                    ...
    ----------------------- ----------------------- ... --------------------------
    2020-01-01 00:00:22.000 2020-01-01 00:00:23.000 ... 45sec. Continuum intensity
    2020-01-01 00:00:22.000 2020-01-01 00:00:23.000 ...         45sec. Magnetogram
    2020-01-01 00:00:22.000 2020-01-01 00:00:23.000 ...         45sec. Dopplergram
    2020-01-01 00:01:07.000 2020-01-01 00:01:08.000 ... 45sec. Continuum intensity
    2020-01-01 00:01:07.000 2020-01-01 00:01:08.000 ...         45sec. Magnetogram
    2020-01-01 00:01:07.000 2020-01-01 00:01:08.000 ...         45sec. Dopplergram
    2020-01-01 00:01:52.000 2020-01-01 00:01:53.000 ... 45sec. Continuum intensity
    2020-01-01 00:01:52.000 2020-01-01 00:01:53.000 ...         45sec. Magnetogram
    2020-01-01 00:01:52.000 2020-01-01 00:01:53.000 ...         45sec. Dopplergram
    2020-01-01 00:02:37.000 2020-01-01 00:02:38.000 ... 45sec. Continuum intensity
    2020-01-01 00:02:37.000 2020-01-01 00:02:38.000 ...         45sec. Magnetogram
    2020-01-01 00:02:37.000 2020-01-01 00:02:38.000 ...         45sec. Dopplergram
    2020-01-01 00:03:22.000 2020-01-01 00:03:23.000 ... 45sec. Continuum intensity
    2020-01-01 00:03:22.000 2020-01-01 00:03:23.000 ...         45sec. Magnetogram
    2020-01-01 00:03:22.000 2020-01-01 00:03:23.000 ...         45sec. Dopplergram
    2020-01-01 00:04:07.000 2020-01-01 00:04:08.000 ... 45sec. Continuum intensity
    2020-01-01 00:04:07.000 2020-01-01 00:04:08.000 ...         45sec. Magnetogram
    2020-01-01 00:04:07.000 2020-01-01 00:04:08.000 ...         45sec. Dopplergram
    2020-01-01 00:04:52.000 2020-01-01 00:04:53.000 ... 45sec. Continuum intensity
    2020-01-01 00:04:52.000 2020-01-01 00:04:53.000 ...         45sec. Magnetogram
    2020-01-01 00:04:52.000 2020-01-01 00:04:53.000 ...         45sec. Dopplergram
    <BLANKLINE>
    <BLANKLINE>

The VSO client returns a lot of information about the records, so the first thing we can do is show only the columns we are interested in.
We can inspect all the available column names in all the responses with the `~.UnifiedResponse.all_colnames` property::

    >>> results.all_colnames  # doctest: +REMOTE_DATA
    ['End Time', 'Extent Length', 'Extent Type', 'Extent Width', 'Info', 'Instrument', 'Physobs', 'Provider', 'Size', 'Source', 'Start Time', 'Wavelength', 'Wavetype', 'fileid']

And then we can pick which ones to see with the :meth:`~.UnifiedResponse.show` method::

    >>> results.show("Start Time", "Instrument", "Physobs", "Wavelength")  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    201 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time       Instrument  Physobs   Wavelength [2]
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
    2020-01-01 00:00:14.000        AIA intensity 1600.0 .. 1600.0
                        ...        ...       ...              ...
    2020-01-01 00:04:47.000        AIA intensity     94.0 .. 94.0
    2020-01-01 00:04:48.000        AIA intensity   335.0 .. 335.0
    2020-01-01 00:04:52.000        AIA intensity 1700.0 .. 1700.0
    2020-01-01 00:04:52.000        AIA intensity   193.0 .. 193.0
    2020-01-01 00:04:53.000        AIA intensity   304.0 .. 304.0
    2020-01-01 00:04:54.000        AIA intensity   131.0 .. 131.0
    2020-01-01 00:04:57.000        AIA intensity   171.0 .. 171.0
    2020-01-01 00:04:57.000        AIA intensity   211.0 .. 211.0
    2020-01-01 00:04:59.000        AIA intensity     94.0 .. 94.0
    2020-01-01 00:05:00.000        AIA intensity   335.0 .. 335.0
    Length = 201 rows
    <BLANKLINE>
    21 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time       Instrument      Physobs        Wavelength [2]
                                                              Angstrom
    ----------------------- ---------- ------------------ ----------------
    2020-01-01 00:00:22.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:00:22.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:00:22.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:01:07.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:01:07.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:01:07.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:01:52.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:01:52.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:01:52.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:02:37.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:02:37.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:02:37.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:03:22.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:03:22.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:03:22.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:04:07.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:04:07.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:04:07.000        HMI       LOS_velocity 6173.0 .. 6174.0
    2020-01-01 00:04:52.000        HMI          intensity 6173.0 .. 6174.0
    2020-01-01 00:04:52.000        HMI LOS_magnetic_field 6173.0 .. 6174.0
    2020-01-01 00:04:52.000        HMI       LOS_velocity 6173.0 .. 6174.0
    <BLANKLINE>
    <BLANKLINE>

To give an example of filtering post-search, let's only return the rows in the table which are line-of-sight magnetograms from HMI or the 94Å passband from AIA.
You can also always do this filtering with the `a.vso.Physobs <sunpy.net.attrs.Physobs>` and `a.Wavelength <sunpy.net.attrs.Wavelength>` attrs in the search command.

First we split the results in to a table for AIA and a table for HMI::

   >>> aia, hmi = results  # doctest: +REMOTE_DATA

We can use boolean indexing to match the value of the ``"Physobs"`` column::

    >>> hmi_los = hmi[hmi["Physobs"] == "LOS_magnetic_field"]  # doctest: +REMOTE_DATA
    >>> hmi_los.show("Start Time", "Instrument", "Wavelength", "Physobs")  # doctest: +REMOTE_DATA
    <sunpy.net.vso.table_response.VSOQueryResponseTable object at ...>
           Start Time       Instrument  Wavelength [2]       Physobs
                                          Angstrom
    ----------------------- ---------- ---------------- ------------------
    2020-01-01 00:00:22.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:01:07.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:01:52.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:02:37.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:03:22.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:04:07.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field
    2020-01-01 00:04:52.000        HMI 6173.0 .. 6174.0 LOS_magnetic_field

To match the ``"Wavelength"`` column we need to account for the fact that VSO results return a wavelength range of ``[min, max]`` so we match the min::

    >>> aia_94 = aia[aia["Wavelength"][:, 0] == 94 * u.AA]  # doctest: +REMOTE_DATA
    >>> aia_94.show("Start Time", "Instrument", "Wavelength", "Physobs")  # doctest: +REMOTE_DATA
    <sunpy.net.vso.table_response.VSOQueryResponseTable object at ...>
           Start Time       Instrument Wavelength [2]  Physobs
                                          Angstrom
    ----------------------- ---------- -------------- ---------
    2020-01-01 00:00:11.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:00:23.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:00:35.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:00:47.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:00:59.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:01:11.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:01:23.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:01:35.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:01:47.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:01:59.000        AIA   94.0 .. 94.0 intensity
                        ...        ...            ...       ...
    2020-01-01 00:03:11.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:03:23.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:03:35.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:03:47.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:03:59.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:04:11.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:04:23.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:04:35.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:04:47.000        AIA   94.0 .. 94.0 intensity
    2020-01-01 00:04:59.000        AIA   94.0 .. 94.0 intensity
    Length = 25 rows

These can then be passed to `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`::

    >>> Fido.fetch(hmi_los, aia_94)  # doctest: +SKIP

.. warning::

   While you can reduce the number of columns and rows in the results, the
   ``fetch()`` method may need certain columns to be present to successfully
   download the files. It is therefore highly recommended that if you are
   planning on downloading data you do not slice out columns, but instead use
   ``.show()`` to only display the ones you are interested in.


.. _downloading_data:

Downloading data
****************
Once you have located your files via a `Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>`, you can download them via `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`::

    >>> downloaded_files = Fido.fetch(results)  # doctest: +SKIP

This downloads the files to the location set in you sunpy config file.
It also returns a `parfive.Results` object ``downloaded_files``, of absolute file paths of where the files have been downloaded to.

You can also specify the path to which you want the data downloaded::

  >>> downloaded_files = Fido.fetch(results, path='/ThisIs/MyPath/to/Data/{file}')  # doctest: +SKIP

This downloads the query results into the directory ``/ThisIs/MyPath/to/Data``, naming each downloaded file with the filename ``{file}`` obtained from the client.
You can also use other properties of the returned query to define the path where the data is saved.
For example, to save the data to a subdirectory named after the instrument, use::

    >>> downloaded_files = Fido.fetch(results, path='./{instrument}/{file}')  # doctest: +SKIP

You can see the list of options that can be specified in path for all the files to be downloaded with ``results.path_format_keys``.

Retrying Downloads
^^^^^^^^^^^^^^^^^^

If any files failed to download, the progress bar will show an incomplete number of files (i.e. 100/150) and the `parfive.Results` object will contain a list of the URLs that failed to transfer and the error associated with them.
This can be accessed with the ``.errors`` attribute or by printing the `~parfive.Results` object::

    >>> print(downloaded_files.errors)  # doctest: +SKIP

The transfer can be retried by passing the `parfive.Results` object back to `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`::

    >>> downloaded_files = Fido.fetch(downloaded_files)  # doctest: +SKIP

doing this will append any newly downloaded file names to the list and replace the ``.errors`` list with any errors that occurred during the second attempt.


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
    CDAWEBClient      Provides access to query and download from the Coordinated Data Analysis Web (CDAWeb).
    EVEClient         Provides access to Level 0CS Extreme ultraviolet Variability Experiment (EVE) data.
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
