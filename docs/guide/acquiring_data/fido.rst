---------------------------------------
Finding and Downloading Data using Fido
---------------------------------------

This guide outlines how to search for and download data using SunPy's
Federated Internet Data Obtainer...or more usually (and
sanely) referred to as Fido.  `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` is a unified interface for searching
and fetching solar physics data irrespective of the underlying
client or webservice through which the data is obtained, e.g. VSO,
JSOC, etc.  It therefore supplies a single, easy and consistent way to
obtain most forms of solar physics data.

Import
------

The `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` object is in
`sunpy.net`. It can be imported as follows::

    >>> from sunpy.net import Fido, attrs as a

Searching for Data Using Fido
-----------------------------

To search for data with Fido, you need to specify attributes to search against.
The (partial) range of allowed attributes are found in the `vso.attrs <sunpy.net.vso.attrs>`
and `jsoc.attrs <sunpy.net.jsoc.attrs>`.
Examples of these attributes are `a.Time <sunpy.net.vso.attrs.Time>`,
`a.Instrument <sunpy.net.vso.attrs.Instrument>`,
`a.Wavelength <sunpy.net.vso.attrs.Wavelength>`, some of these attributes are
client specific, such as `a.vso <sunpy.net.vso.attrs>` or
`a.jsoc <sunpy.net.jsoc.attrs>`.::

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('lyra'))

This returns an `~sunpy.net.fido_factory.UnifiedResponse` object containing
information on the available online files which fit the criteria specified by
the attrs objects in the above call. It does not download the files. For
instructions on how to download data using Fido, see :ref:`downloading_data`.

To see a summary of results of our query, simple type the name of the
variable set to the Fido search, in this case, result::

    >>> result  # doctest: +SKIP
    <sunpy.net.fido_factory.UnifiedResponse object at 0x7fe6258ab630>
    Results from 1 Provider:

    3 Results from the LYRAClient:
        Start Time           End Time      Source Instrument Wavelength
          str19               str19         str6     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2012-03-04 00:00:00 2012-03-06 00:00:00 Proba2       lyra        nan
    2012-03-04 00:00:00 2012-03-06 00:00:00 Proba2       lyra        nan
    2012-03-04 00:00:00 2012-03-06 00:00:00 Proba2       lyra        nan

Queries can be made more flexible or specific by adding more attrs objects to
the `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` search. Specific
passbands can be searched for by supplying an `~astropy.units.Quantity` to the
`a.Wavelength <sunpy.net.vso.attrs.Wavelength>` attribute::

    >>> import astropy.units as u
    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('norh'),
    ...                      a.Wavelength(17*u.GHz))  # doctest: +REMOTE_DATA

Data of a given cadence can also be specified using the Sample attribute. To
search for data at a given cadence use the
`a.vso.Sample <sunpy.net.vso.attrs.Sample>` attribute.
`a.vso.Sample <sunpy.net.vso.attrs.Sample>` is only supported by the
`sunpy.net.vso.VSOClient` hence it has the ``a.vso`` prefix. Attributes
like this which are client specific will result in
`Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` only searching that
client for results, in this case VSO.::

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('aia'),
    ...                      a.Wavelength(171*u.angstrom), a.vso.Sample(10*u.minute))  # doctest: +REMOTE_DATA

To search for data from multiple instruments, wavelengths, times etc., use the
pipe ``|`` operator. This joins queries together just as the logical ``OR``
operator would::

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'),
    ...                      a.Instrument('lyra') | a.Instrument('rhessi'))  # doctest: +REMOTE_DATA

    >>> result = Fido.search(a.Time('2012/3/4', '2012/3/6'), a.Instrument('aia'),
    ...                      a.Wavelength(171*u.angstrom) | a.Wavelength(94*u.angstrom))  # doctest: +REMOTE_DATA


Indexing search results
-----------------------

The `~sunpy.net.fido_factory.UnifiedResponse` that Fido returns can be
indexed to access a subset of the search results. When doing this, the
results should be treated as a two-dimensional array in which the first
dimension corresponds to the clients which have returned results and the
second to the records returned.

For example, the following code returns a response containing LYRA data from
the `~sunpy.net.dataretriever.sources.LYRAClient`, and EVE data from the
`~sunpy.net.vso.vso.VSOClient`::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"),
    ...                       a.Instrument("lyra") | a.Instrument("eve"))  # doctest: +REMOTE_DATA

If you then wanted to inspect just the LYRA data for the whole time range
specified in the search, you would index this response to see just the
results returned by the `~sunpy.net.dataretriever.sources.LYRAClient`::

    >>> results[0, :]  # doctest: +SKIP
    <sunpy.net.fido_factory.UnifiedResponse object at 0x7fe61fdf1b00>
    Results from 1 Provider:

    2 Results from the LYRAClient:
        Start Time           End Time      Source Instrument Wavelength
          str19               str19         str6     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2012-01-01 00:00:00 2012-01-02 00:00:00 Proba2       lyra        nan
    2012-01-01 00:00:00 2012-01-02 00:00:00 Proba2       lyra        nan


Or, equivalently::

    >>> results[0]  # doctest: +SKIP
    <sunpy.net.fido_factory.UnifiedResponse object at 0x7fe625811748>
    Results from 1 Provider:

    2 Results from the LYRAClient:
        Start Time           End Time      Source Instrument Wavelength
          str19               str19         str6     str4       str3
    ------------------- ------------------- ------ ---------- ----------
    2012-01-01 00:00:00 2012-01-02 00:00:00 Proba2       lyra        nan
    2012-01-01 00:00:00 2012-01-02 00:00:00 Proba2       lyra        nan


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
----------------
Once you have located your files via a
`Fido.search <sunpy.net.fido_factory.UnifiedDownloaderFactory.search>`, you can
download them via `Fido.fetch <sunpy.net.fido_factory.UnifiedDownloaderFactory.fetch>`::

    >>> downloaded_files = Fido.fetch(results)  # doctest: +SKIP

This downloads the files to the location set in you sunpy config
file.  It also returns a list ``downloaded_files``, of absolute file paths
of where the files have been downloaded to.

You can also specify the path to which you want the data downloaded::

  >>> downloaded_files = Fido.fetch(results, path='/ThisIs/MyPath/to/Data/{file}.fits')  # doctest: +SKIP

This downloads the query results into the directory
``/ThisIs/MyPath/to/Data``, naming each downloaded file with the
filename ``{file}`` obtained from the client, and appended with the suffix
``.fits``. You can also use other properties of the returned query
to define the path where the data is saved.  For example, to save the
data to a subdirectory named after the instrument, use

    >>> downloaded_files = Fido.fetch(results, path='./{instrument}/{file}.fits')  # doctest: +SKIP

You can see the list of options that can be specified in path for all the files
to be downloaded with ``results.response_block_properties``.
