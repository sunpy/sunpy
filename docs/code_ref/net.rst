SunPy net
=========

SunPy's net submodule contains a lot of different code for accessing various
solar physics related web services. This submodule contains many layers. Most
users should use `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`, which
is an interface to multiple sources including all the sources implemented in
`~sunpy.net.dataretriever` as well as `~sunpy.net.vso` and `~sunpy.net.jsoc`.
`Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>` can be used like so::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Instrument('lyra'))  # doctest: +REMOTE_DATA
        [<class 'sunpy.net.dataretriever.client.QueryResponse'><Table length=2>
             Start Time           End Time      Source Instrument Wavelength
               str19               str19         str6     str4       str3
        ------------------- ------------------- ------ ---------- ----------
        2012-01-01 00:00:00 2012-01-02 00:00:00 Proba2       lyra        nan
        2012-01-01 00:00:00 2012-01-02 00:00:00 Proba2       lyra        nan]
    >>> files = Fido.fetch(results)  # doctest: +SKIP

.. automodapi:: sunpy.net
   :no-heading:

.. automodapi:: sunpy.net.fido_factory

Dataretriever
-------------

.. automodapi:: sunpy.net.dataretriever
   :allowed-package-names: sources
   :headings: ^#

.. automodapi:: sunpy.net.dataretriever.sources
   :headings: #~


VSO
---

.. automodapi:: sunpy.net.vso
   :headings: ^#

.. automodapi:: sunpy.net.vso.attrs
   :headings: #~


HEK
---

.. automodapi:: sunpy.net.hek
    :headings: ^#

.. automodapi:: sunpy.net.hek2vso
    :headings: ^#


HELIO
-----

.. automodapi:: sunpy.net.helio
    :headings: ^#

.. automodapi:: sunpy.net.helio.hec
    :headings: #~


JSOC
----

.. automodapi:: sunpy.net.jsoc
    :headings: ^#

.. automodapi:: sunpy.net.jsoc.attrs
    :headings: #~
