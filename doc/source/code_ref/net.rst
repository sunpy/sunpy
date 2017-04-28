SunPy net
=========

SunPy's net submodule contains a lot of different code for accessing various
solar physics related web services. This submodule contains many layers. Most
users should use ``Fido``, which is an interface to multiple sources including
all the sources implemented in `~sunpy.net.dataretriever` as well as
`~sunpy.net.vso` and `~sunpy.net.jsoc`. ``Fido`` can be used like so::

>>> from sunpy.net import Fido, attrs as a
>>> results = Fido.search(a.Time("2012/1/1", "2012/1/2"), a.Instrument('lyra'))
>>> files = Fido.fetch(results)

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

