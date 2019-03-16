SunPy net
=========

.. toctree::
   :maxdepth: 1

   fido
   hek
   jsoc
   dataretriever
   attr

SunPy's net submodule contains a lot of different code for accessing various
solar physics related web services. This submodule contains many layers. Most
users should use `Fido <sunpy.net.fido_factory.UnifiedDownloaderFactory>`, which
is an interface to multiple sources including all the sources implemented in
`~sunpy.net.dataretriever` as well as `~sunpy.net.vso` and `~sunpy.net.jsoc`.
Fido ~`sunpy.net.fido_factory.UnifiedDownloaderFactory` can be used like so::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Instrument('lyra'))  # doctest: +REMOTE_DATA
    >>> files = Fido.fetch(results)  # doctest: +SKIP

.. automodapi:: sunpy.net
   :no-heading:
