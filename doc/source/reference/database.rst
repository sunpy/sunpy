Database
========

..
    .. currentmodule:: sunpy.database
    TODO: - WHY was this package developed? WHAT are the use-cases?
    - WHICH methods can be undone?

.. automodule:: sunpy.database

Module documentation
^^^^^^^^^^^^^^^^^^^^

The Database class
""""""""""""""""""

.. autoclass:: sunpy.database.Database
   :members:

Exceptions
""""""""""
.. class:: sunpy.database.EntryAlreadyAddedError

    This exception is raised if a database entry is attempted to be added
    to the database although it was already saved in it.

.. class:: sunpy.database.EntryAlreadyStarredError

    This exception is raised if a database entry is marked as starred
    using :meth:`Database.star` although it was already starred before
    this operation.

.. class:: sunpy.database.EntryAlreadyUnstarredError

    This exception is raised if the star mark from a database entry is
    attempted to be removed although the entry is not starred.

.. class:: sunpy.database.EntryNotFoundError

    This exception is raised if a database entry cannot be found by its
    unique ID.

.. class:: sunpy.database.TagAlreadyAssignedError

    This exception is raised if it is attempted to assign a tag to a
    database entry but the database entry already has this tag assigned.

.. class:: sunpy.database.NoSuchTagError

    This exception is raised if a tag cannot be found in a database by its
    name.

.. FIXME: to be documented
    .. class:: sunpy.database.commands.NoSuchEntryError
    .. class:: sunpy.database.commands.EmptyCommandStackError
..
    tables
    """"""
    .. autoclass:: sunpy.database.tables.DatabaseEntry
       :members:
    .. autoclass:: sunpy.database.tables.FitsHeaderEntry
       :members:
    .. autoclass:: sunpy.database.tables.Tag
       :members:
    utility functions
    """""""""""""""""
    .. autofunction:: sunpy.database.tables.entries_from_query_result
    .. autofunction:: sunpy.database.tables.entries_from_path
    Caching
    """""""
    .. autoclass:: sunpy.database.caching.BaseCache
       :members: callback
    .. autoclass:: sunpy.database.caching.LRUCache
    .. autoclass:: sunpy.database.caching.LFUCache
