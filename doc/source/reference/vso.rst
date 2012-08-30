.. _vso:

===
VSO
===

.. currentmodule:: sunpy.net.vso

Introduction
^^^^^^^^^^^^
The VSO can be queried in two different ways: the legacy or the standard
API. The legacy API tries to mimic the IDL vso_search as closely as
possible (and sane for Python), the standard API takes a new approach
and allows complex queries to be easily created.

There also exist two different clients for the VSO. The `VSOClient` is designed 
for computer programs while `InteractiveVSOClient` is designed for use in the 
REPL. `InteractiveVSOClient` queries the user to provide missing information 
using the keyboard while `VSOClient` can be programmed to supply missing 
information by subclassing and overriding methods.

`InteractiveVSOClient` also offeres a search method that can be used to do both 
legacy and standard queries with one method.

A global `search` and `get` method is exposed on module level that operate on 
a default `InteractiveVSOClient`.

Standard Queries
^^^^^^^^^^^^^^^^
`VSOClient.query` takes an arbitrary amount of Attrs which are automatically 
ANDed together for convenience, `a | b` can be used for OR::

   client = vso.VSOClient()
   result = client.query(
      vso.attrs.Time((2010, 1, 1), (2010, 1, 1, 1)),
      vso.attrs.Instrument('eit') | vso.attrs.Instrument('sot')
   )

Available attributes
--------------------
.. automodule:: sunpy.net.vso.attrs
   :members:
   :undoc-members:
   :exclude-members: Range

Legacy Queries
^^^^^^^^^^^^^^
See `VSOClient.query_legacy`

Module documentation
^^^^^^^^^^^^^^^^^^^^
.. automodule:: sunpy.net.vso
   :members: