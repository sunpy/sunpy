===
VSO
===

Introduction
^^^^^^^^^^^^
The VSO can be queried in two different ways: the legacy or the standard
API. The legacy API tries to mimick the IDL vso_search as closely as
possible (and sane for Python), the standard API takes a new approach
and allows complex queries to be easily created.

There also exist two different clients for the VSO. The :py:class:`VSOClient`
is designed for computer programs while :py:class:`InteractiveVSOClient` is
designed for use in the REPL. :py:class:`InteractiveVSOClient` queries the
user to provide missing information using the keyboard while
:py:class:`VSOClient` can be programmed to supply missing information by
subclassing and overriding methods.

:py:class:`InteractiveVSOClient` also offeres a search method that can
be used to do both legacy and standard queries with one method.

Standard Queries
^^^^^^^^^^^^^^^^
:py:meth:`VSOClient.query` takes an arbitrary amount of Attrs which are
automatically ANDed together for convenience::

   client = vso.VSOClient()
   client.query(vso.Time.dt((2010, 1, 1), (2010, 1, 1, 1)), vso.Instrument('eit') | vso.Instrument('ait'))

Legacy Queries
^^^^^^^^^^^^^^

Module documentation
^^^^^^^^^^^^^^^^^^^^
.. automodule:: sunpy.net.vso
   :members:
