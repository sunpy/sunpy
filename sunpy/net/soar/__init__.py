"""
A Fido client for accessing data in the `Solar Orbiter Archive (SOAR) <https://soar.esac.esa.int/soar/>`__.

.. note::

    This client used to live in the ``sunpy_soar`` package, which is no longer needed, you should uninstall or not import ``sunpy_soar`` when running sunpy >= 8.0.

Examples
--------

>>> import sunpy.net.attrs as a
>>> from sunpy.net import Fido

>>> instrument = a.Instrument("EUI")
>>> time = a.Time("2021-02-01", "2021-02-01 01:00")
>>> level = a.Level(1)
>>> product = a.soar.Product("EUI-FSI174-IMAGE")

>>> Fido.search(instrument & time & level & product) # doctest: +REMOTE_DATA
<sunpy.net.fido_factory.UnifiedResponse object at ...>
Results from 1 Provider:
<BLANKLINE>
1 Results from the SOARClient:
<BLANKLINE>
Instrument   Data product   Level        Start time               End time        Filesize SOOP Name Detector Wavelength
                                                                                   Mbyte
---------- ---------------- ----- ----------------------- ----------------------- -------- --------- -------- ----------
       EUI eui-fsi174-image    L1 2021-02-01 00:45:12.228 2021-02-01 00:45:22.228    3.393      none      FSI      174.0
<BLANKLINE>
<BLANKLINE>

"""

from sunpy.net.soar.client import SOARClient
