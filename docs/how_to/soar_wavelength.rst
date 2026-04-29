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

    >>> instrument = a.Instrument("METIS")
    >>> time = a.Time("2023-02-01 01:00", "2023-02-01 05:00")
    >>> level = a.Level(2)
    >>> wavelength = a.Wavelength(121.6 * u.AA)
    >>> result = Fido.search(instrument & time & level & wavelength) # doctest: +REMOTE_DATA
    >>> result  # doctest: +REMOTE_DATA
    <sunpy.net.fido_factory.UnifiedResponse object at ...>
    Results from 1 Provider:
    <BLANKLINE>
    ... Results from the SOARClient:
    <BLANKLINE>
    Instrument  Data product  Level        Start time               End time        Filesize SOOP Name Detector Wavelength
                                                                                     Mbyte
    ---------- -------------- ----- ----------------------- ----------------------- -------- --------- -------- ----------
         METIS metis-uv-image    L2 2023-02-01 01:00:48.690 2023-02-01 01:11:46.866     0.85      none      UVD      121.6
    ...


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
    ... Results from the SOARClient:
    <BLANKLINE>
    Instrument    Data product    Level        Start time               End time        Filesize SOOP Name Detector Wavelength
                                                                                         Mbyte
    ---------- ------------------ ----- ----------------------- ----------------------- -------- --------- -------- ----------
         METIS        metis-vl-tb    L2 2023-02-01 01:00:00.147 2023-02-01 01:24:37.923    12.64      none      VLD      610.0
    ...
