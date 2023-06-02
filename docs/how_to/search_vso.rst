.. _sunpy-how-to-search-the-vso:

******************************************
Search the Virtual Solar Observatory (VSO)
******************************************

To search the Virtual Solar Observatory (VSO) for SDO AIA data in all channels over a given time range, use the timerange (`~sunpy.net.attrs.Time`) and the instrument (`~sunpy.net.attrs.Instrument`) attrs,

.. code-block:: python

    >>> import astropy.units as u

    >>> from sunpy.net import Fido
    >>> from sunpy.net import attrs as a

    >>> time_range = a.Time('2020/03/04 00:00:15', '2020/03/04 00:00:30')
    >>> result = Fido.search(time_range, a.Instrument.aia)  # doctest: +REMOTE_DATA
    >>> result  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    11 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 745.677 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source  ... Extent Type   Size
                                                            ...              Mibyte
    ----------------------- ----------------------- ------  ... ----------- --------
    2020-03-04 00:00:16.000 2020-03-04 00:00:17.000    SDO  ...    FULLDISK 64.64844
    2020-03-04 00:00:17.000 2020-03-04 00:00:18.000    SDO  ...    FULLDISK 64.64844
    2020-03-04 00:00:18.000 2020-03-04 00:00:19.000    SDO  ...    FULLDISK 64.64844
                        ...                     ...    ...  ...         ...      ...
    2020-03-04 00:00:28.000 2020-03-04 00:00:29.000    SDO  ...    FULLDISK 64.64844
    2020-03-04 00:00:29.000 2020-03-04 00:00:30.000    SDO  ...    FULLDISK 64.64844
    2020-03-04 00:00:30.000 2020-03-04 00:00:31.000    SDO  ...    FULLDISK 64.64844
    <BLANKLINE>
    <BLANKLINE>

You can also search for specific physical observables using `~sunpy.net.attrs.Physobs`.
For example, you can search for line-of-sight (LOS) magnetic field measurements from HMI in the same time range,

.. code-block:: python

    >>> Fido.search(time_range, a.Instrument.hmi, a.Physobs.los_magnetic_field)  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    1 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2020-03-04 00:00:26.000 2020-03-04 00:00:27.000    SDO ...    FULLDISK -0.00098
    <BLANKLINE>
    <BLANKLINE>

You can also use relational operators when constructing queries.
For example, the AIA query above can also be expressed using the AND (&) operator,

.. code-block:: python

    >>> Fido.search(time_range & a.Instrument.aia)  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    11 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 745.677 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2020-03-04 00:00:16.000 2020-03-04 00:00:17.000    SDO ...    FULLDISK 64.64844
    2020-03-04 00:00:17.000 2020-03-04 00:00:18.000    SDO ...    FULLDISK 64.64844
    2020-03-04 00:00:18.000 2020-03-04 00:00:19.000    SDO ...    FULLDISK 64.64844
    2020-03-04 00:00:21.000 2020-03-04 00:00:22.000    SDO ...    FULLDISK 64.64844
    2020-03-04 00:00:21.000 2020-03-04 00:00:22.000    SDO ...    FULLDISK 64.64844
    2020-03-04 00:00:23.000 2020-03-04 00:00:24.000    SDO ...    FULLDISK 64.64844
    2020-03-04 00:00:24.000 2020-03-04 00:00:25.000    SDO ...    FULLDISK 64.64844
    2020-03-04 00:00:28.000 2020-03-04 00:00:29.000    SDO ...    FULLDISK 64.64844
    2020-03-04 00:00:28.000 2020-03-04 00:00:29.000    SDO ...    FULLDISK 64.64844
    2020-03-04 00:00:29.000 2020-03-04 00:00:30.000    SDO ...    FULLDISK 64.64844
    2020-03-04 00:00:30.000 2020-03-04 00:00:31.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    <BLANKLINE>

Additionally, multiple operators can be chained together to, for example, search for only the 171 channel,

.. code-block:: python

    >>> Fido.search(time_range & a.Instrument.aia & a.Wavelength(171*u.angstrom))  # doctest: +REMOTE_DATA
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
    2020-03-04 00:00:21.000 2020-03-04 00:00:22.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    <BLANKLINE>

The OR operator (``|``) can also be used to construct queries.
For example, to search for AIA data in this same time range from both the 94 and 171 channels,

.. code-block:: python

    >>> Fido.search(time_range,
    ...             a.Instrument.aia,
    ...             a.Wavelength(171*u.angstrom) | a.Wavelength(94*u.angstrom))  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 2 Providers:
    <BLANKLINE>
    1 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 67.789 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2020-03-04 00:00:21.000 2020-03-04 00:00:22.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    1 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 67.789 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2020-03-04 00:00:23.000 2020-03-04 00:00:24.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    <BLANKLINE>

These relational operators are particularly useful when searching a given time interval for multiple instruments.
For example, to find the HMI LOS magnetic field data and the AIA 94 and 171 data in the given time interval,

.. code-block:: python

    >>> aia_params = a.Instrument.aia & (a.Wavelength(171*u.angstrom) | a.Wavelength(94*u.angstrom))
    >>> hmi_params = a.Instrument.hmi & a.Physobs.los_magnetic_field
    >>> Fido.search(time_range, aia_params | hmi_params)  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 3 Providers:
    <BLANKLINE>
    1 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 67.789 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2020-03-04 00:00:21.000 2020-03-04 00:00:22.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    1 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    Total estimated size: 67.789 Mbyte
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2020-03-04 00:00:23.000 2020-03-04 00:00:24.000    SDO ...    FULLDISK 64.64844
    <BLANKLINE>
    1 Results from the VSOClient:
    Source: http://vso.stanford.edu/cgi-bin/search
    <BLANKLINE>
           Start Time               End Time        Source ... Extent Type   Size
                                                           ...              Mibyte
    ----------------------- ----------------------- ------ ... ----------- --------
    2020-03-04 00:00:26.000 2020-03-04 00:00:27.000    SDO ...    FULLDISK -0.00098
    <BLANKLINE>
    <BLANKLINE>
