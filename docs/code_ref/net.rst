*************************
Remote data (`sunpy.net`)
*************************

`sunpy.net` contains a lot of different classess for accessing various solar physics related web services.

Most users should use `~sunpy.net.Fido`, which is an interface to multiple sources including all the sources implemented in `~sunpy.net.dataretriever` as well as `~sunpy.net.vso` and `~sunpy.net.jsoc`.

.. automodapi:: sunpy.net
   :include-all-objects:

.. automodapi:: sunpy.net.fido_factory

.. automodapi:: sunpy.net.attrs

.. automodapi:: sunpy.net.vso

.. automodapi:: sunpy.net.vso.attrs

.. automodapi:: sunpy.net.jsoc

.. automodapi:: sunpy.net.jsoc.attrs

.. automodapi:: sunpy.net.hek

.. automodapi:: sunpy.net.hek.attrs

.. automodapi:: sunpy.net.hek2vso

.. automodapi:: sunpy.net.cdaweb

.. automodapi:: sunpy.net.helio

.. automodapi:: sunpy.net.helio.attrs

.. automodapi:: sunpy.net.dataretriever

.. automodapi:: sunpy.net.dataretriever.attrs

.. automodapi:: sunpy.net.helioviewer

Internal Classes and Functions
==============================

These classes and functions are designed to be used to help develop new clients for `sunpy.net.Fido`.

.. automodapi:: sunpy.net.base_client

.. automodapi:: sunpy.net.dataretriever.client

.. automodapi:: sunpy.net.attr

.. automodapi:: sunpy.net.scraper
