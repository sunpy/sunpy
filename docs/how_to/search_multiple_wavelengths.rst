.. _sunpy-how-to-search-for-multiple-wavelengths-with-fido:

*****************************************
Search for multiple wavelengths with Fido
*****************************************

Use the `~sunpy.net.attrs.Wavelength` to search for a particular wavelength:

.. code-block:: python

    >>> from astropy import units as u

    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a

    >>> time_range = a.Time("2022-02-20 00:00:00", "2022-02-20 00:00:30")
    >>> aia_search = Fido.search(time_range,
    ...                          a.Instrument.aia,
    ...                          a.Wavelength(171*u.angstrom))  # doctest: +REMOTE_DATA
    >>> aia_search  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    2 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 135.578 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2022-02-20 00:00:09.000 2022-02-20 00:00:10.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:21.000 2022-02-20 00:00:22.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    <BLANKLINE>

The "|" operator can be used to combine multiple wavelengths:

.. code-block:: python

    >>> aia_search = Fido.search(time_range,
    ...                          a.Instrument.aia,
    ...                          a.Wavelength(171*u.angstrom) | a.Wavelength(193*u.angstrom))  # doctest: +REMOTE_DATA
    >>> aia_search  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    2 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 135.578 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2022-02-20 00:00:09.000 2022-02-20 00:00:10.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:21.000 2022-02-20 00:00:22.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    3 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 203.366 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2022-02-20 00:00:04.000 2022-02-20 00:00:05.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:16.000 2022-02-20 00:00:17.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:28.000 2022-02-20 00:00:29.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    <BLANKLINE>

When searching for more than two wavelengths, it is more practical to use the :func:`sunpy.net.attrs.AttrOr` function:

.. code-block:: python

    >>> wavelengths = [94, 131, 171, 193, 211]*u.angstrom
    >>> aia_search = Fido.search(time_range,
    ...                         a.Instrument.aia,
    ...                         a.AttrOr([a.Wavelength(wav) for wav in wavelengths]))  # doctest: +REMOTE_DATA
    >>> aia_search  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 5 Providers:
    <BLANKLINE>
    2 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 135.578 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2022-02-20 00:00:11.000 2022-02-20 00:00:12.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:23.000 2022-02-20 00:00:24.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    3 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 203.366 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2022-02-20 00:00:06.000 2022-02-20 00:00:07.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:18.000 2022-02-20 00:00:19.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:30.000 2022-02-20 00:00:31.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    2 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 135.578 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2022-02-20 00:00:09.000 2022-02-20 00:00:10.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:21.000 2022-02-20 00:00:22.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    3 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 203.366 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2022-02-20 00:00:04.000 2022-02-20 00:00:05.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:16.000 2022-02-20 00:00:17.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:28.000 2022-02-20 00:00:29.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    2 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 135.578 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2022-02-20 00:00:09.000 2022-02-20 00:00:10.000    SDO ...    FULLDISK 64.64844
    2022-02-20 00:00:21.000 2022-02-20 00:00:22.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    <BLANKLINE>
