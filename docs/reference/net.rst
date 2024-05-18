Remote data (`sunpy.net`)
*************************

``sunpy.net`` contains a lot of different code for accessing various solar
physics related web services. This submodule contains many layers. Most users
should use `~sunpy.net.Fido`, which is an interface to multiple sources
including all the sources implemented in `~sunpy.net.dataretriever` as well as
`~sunpy.net.vso` and `~sunpy.net.jsoc`. `~sunpy.net.Fido` can be used like so::

    >>> from sunpy.net import Fido, attrs as a
    >>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Instrument.lyra)  # doctest: +REMOTE_DATA
    >>> files = Fido.fetch(results)  # doctest: +SKIP

.. automodapi:: sunpy.net
   :include-all-objects:
   :headings: =-

.. automodapi:: sunpy.net.attrs
   :headings: ^"

.. automodapi:: sunpy.net.fido_factory
   :headings: ^"

VSO
---

.. automodapi:: sunpy.net.vso
   :headings: ^"

.. automodapi:: sunpy.net.vso.attrs
   :headings: ^"

Dataretriever
-------------

.. automodapi:: sunpy.net.dataretriever
   :headings: ^"

.. automodapi:: sunpy.net.dataretriever.attrs.goes
   :headings: ^"

JSOC
----

.. automodapi:: sunpy.net.jsoc
   :headings: ^"

.. automodapi:: sunpy.net.jsoc.attrs
   :headings: ^"

HEK
---

.. automodapi:: sunpy.net.hek
   :headings: ^"

.. automodapi:: sunpy.net.hek.attrs
   :headings: ^"

.. automodapi:: sunpy.net.hek2vso
   :headings: ^"

CDAWeb
------

.. automodapi:: sunpy.net.cdaweb
   :headings: ^"


HELIO
-----

.. automodapi:: sunpy.net.helio
   :headings: ^"

.. automodapi:: sunpy.net.helio.attrs
   :headings: ^"


Internal Classes and Functions
==============================

These classes and functions are designed to be used to help develop new clients
for `sunpy.net.Fido`.

.. automodapi:: sunpy.net.base_client

.. automodapi:: sunpy.net.dataretriever.client

.. automodapi:: sunpy.net.attr

.. automodapi:: sunpy.net.scraper

.. automodapi:: sunpy.net.scraper_utils