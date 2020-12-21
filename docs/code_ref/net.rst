*********
SunPy net
*********

SunPy's net submodule contains a lot of different code for accessing various
solar physics related web services. This submodule contains many layers. Most
users should use `~sunpy.net.Fido`, which
is an interface to multiple sources including all the sources implemented in
`~sunpy.net.dataretriever` as well as `~sunpy.net.vso` and `~sunpy.net.jsoc`.
`~sunpy.net.Fido` can be used like so::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Instrument.lyra)  # doctest: +REMOTE_DATA
    >>> files = Fido.fetch(results)  # doctest: +SKIP

.. automodapi:: sunpy.net
   :include-all-objects:

.. automodapi:: sunpy.net.attrs

.. automodapi:: sunpy.net.fido_factory


VSO
===

.. automodapi:: sunpy.net.vso

.. automodapi:: sunpy.net.vso.attrs


Dataretriever
=============

.. automodapi:: sunpy.net.dataretriever

.. automodapi:: sunpy.net.dataretriever.attrs.goes

JSOC
====

.. automodapi:: sunpy.net.jsoc

.. automodapi:: sunpy.net.jsoc.attrs


HEK
===

.. automodapi:: sunpy.net.hek

.. automodapi:: sunpy.net.hek2vso


HELIO
=====

.. automodapi:: sunpy.net.helio

.. automodapi:: sunpy.net.helio.hec

Helioviewer
===========

.. automodapi:: sunpy.net.helioviewer


Internal Classes and Functions
==============================

These classes and functions are designed to be used to help develop new clients
for `sunpy.net.Fido`.

.. automodapi:: sunpy.net.base_client

.. automodapi:: sunpy.net.attr
