.. _sunpy-soar-how-to-query-wavelength:

***************************************************************
Query for Solar Orbiter data using the ``Wavelength`` attribute
***************************************************************

``sunpy-soar`` provides convenient methods to construct queries using ``a.Wavelength`` for several different remote sensing instruments available through the SOAR.
In this guide, we will demonstrate how we can query data using ``a.Wavelength`` for different instruments.
For instruments EUI, METIS and SOLOHI we get query results in form of wavelength.
However, at this time we cannot search for wavelengths for instruments SPICE, PHI, and STIX as this information is not yet available in the SOAR.

For instruments EUI, METIS and SOLOHI passing a single Wavelength
=================================================================

When a single wavelength is provided it is interpreted as the wavelength.

.. code-block:: python

    >>> import astropy.units as u
    >>> import sunpy.net.attrs as a
    >>> from sunpy.net import Fido
    >>> import sunpy.net.soar

    >>> instrument = a.Instrument("METIS")
    >>> time = a.Time("2023-02-01 01:00", "2023-02-01 05:00")
    >>> level = a.Level(2)
    >>> wavelength = a.Wavelength(121.6 * u.AA)
    >>> result = Fido.search(instrument & time & level & wavelength) # doctest: +REMOTE_DATA
    >>> result  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    8 Results from the SOARClient:
    <BLANKLINE>
    Instrument  Data product  Level        Start time               End time        Filesize SOOP Name Detector Wavelength
                                                                                     Mbyte
    ---------- -------------- ----- ----------------------- ----------------------- -------- --------- -------- ----------
         METIS metis-uv-image    L2 2023-02-01 01:00:48.690 2023-02-01 01:11:46.866     0.85      none      UVD      121.6
         METIS metis-uv-image    L2 2023-02-01 01:30:48.680 2023-02-01 01:41:45.540     0.85      none      UVD      121.6
         METIS metis-uv-image    L2 2023-02-01 02:00:48.671 2023-02-01 02:11:44.213    12.64      none      UVD      121.6
         METIS metis-uv-image    L2 2023-02-01 02:30:48.661 2023-02-01 02:41:42.887    12.64      none      UVD      121.6
         METIS metis-uv-image    L2 2023-02-01 03:00:48.652 2023-02-01 03:11:41.560    12.64      none      UVD      121.6
         METIS metis-uv-image    L2 2023-02-01 03:30:48.643 2023-02-01 03:41:40.233    12.64      none      UVD      121.6
         METIS metis-uv-image    L2 2023-02-01 04:00:48.633 2023-02-01 04:11:38.907    12.64      none      UVD      121.6
         METIS metis-uv-image    L2 2023-02-01 04:30:38.163 2023-02-01 04:40:37.625    12.64      none      UVD      121.6
    <BLANKLINE>
    <BLANKLINE>

For instruments EUI, METIS and SOLOHI passing a range of Wavelength
===================================================================

When a range of wavelength is provided, it is interpreted as the wavemin and wavemax.

.. code-block:: python

    >>> wavelength = a.Wavelength(580 * u.AA, 640 * u.AA)
    >>> result = Fido.search(instrument & time & level & wavelength) # doctest: +REMOTE_DATA
    >>> result  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    64 Results from the SOARClient:
    <BLANKLINE>
    Instrument    Data product    Level        Start time               End time        Filesize SOOP Name Detector Wavelength
                                                                                         Mbyte
    ---------- ------------------ ----- ----------------------- ----------------------- -------- --------- -------- ----------
         METIS        metis-vl-tb    L2 2023-02-01 01:00:00.147 2023-02-01 01:24:37.923    12.64      none      VLD      610.0
         METIS    metis-vl-stokes    L2 2023-02-01 01:00:00.147 2023-02-01 01:24:37.923   21.067      none      VLD      610.0
         METIS        metis-vl-pb    L2 2023-02-01 01:00:00.147 2023-02-01 01:24:37.923    12.64      none      VLD      610.0
         METIS     metis-vl-image    L2 2023-02-01 01:00:00.147 2023-02-01 01:23:05.525   12.643      none      VLD      610.0
         METIS metis-vl-pol-angle    L2 2023-02-01 01:00:00.147 2023-02-01 01:24:37.923    12.64      none      VLD      610.0
         METIS     metis-vl-image    L2 2023-02-01 01:00:30.944 2023-02-01 01:23:36.322   12.643      none      VLD      610.0
         METIS     metis-vl-image    L2 2023-02-01 01:01:01.741 2023-02-01 01:24:07.127   12.643      none      VLD      610.0
         METIS     metis-vl-image    L2 2023-02-01 01:01:32.538 2023-02-01 01:24:37.922   12.643      none      VLD      610.0
         METIS    metis-vl-stokes    L2 2023-02-01 01:30:00.150 2023-02-01 01:54:37.922   21.067      none      VLD      610.0
         METIS metis-vl-pol-angle    L2 2023-02-01 01:30:00.150 2023-02-01 01:54:37.922    12.64      none      VLD      610.0
           ...                ...   ...                     ...                     ...      ...       ...      ...        ...
         METIS     metis-vl-image    L2 2023-02-01 04:01:01.774 2023-02-01 04:24:07.158   50.388      none      VLD      610.0
         METIS     metis-vl-image    L2 2023-02-01 04:01:32.571 2023-02-01 04:24:37.955   50.388      none      VLD      610.0
         METIS metis-vl-pol-angle    L2 2023-02-01 04:30:00.201 2023-02-01 04:54:37.981   50.388      none      VLD      610.0
         METIS    metis-vl-stokes    L2 2023-02-01 04:30:00.201 2023-02-01 04:54:37.981   83.981      none      VLD      610.0
         METIS        metis-vl-tb    L2 2023-02-01 04:30:00.201 2023-02-01 04:54:37.981   50.388      none      VLD      610.0
         METIS        metis-vl-pb    L2 2023-02-01 04:30:00.201 2023-02-01 04:54:37.981   50.388      none      VLD      610.0
         METIS     metis-vl-image    L2 2023-02-01 04:30:00.201 2023-02-01 04:53:05.585   50.388      none      VLD      610.0
         METIS     metis-vl-image    L2 2023-02-01 04:30:30.999 2023-02-01 04:53:36.383   50.388      none      VLD      610.0
         METIS     metis-vl-image    L2 2023-02-01 04:31:01.796 2023-02-01 04:54:07.184   50.388      none      VLD      610.0
         METIS     metis-vl-image    L2 2023-02-01 04:31:32.593 2023-02-01 04:54:37.979   50.388      none      VLD      610.0
    Length = 64 rows
    <BLANKLINE>
    <BLANKLINE>
