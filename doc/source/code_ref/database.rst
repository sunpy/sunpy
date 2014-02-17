SunPy Database
==============
..
    TODO:
        - WHY was this package developed?
        - WHAT are the use-cases?
        - WHICH methods can be undone?
        - WHICH methods are the most important ones to know?

.. automodapi:: sunpy.database
    :no-inheritance-diagram:

Exceptions
^^^^^^^^^^
.. class:: sunpy.database.EntryAlreadyAddedError(database_entry)

    This exception is raised if a database entry is attempted to be added
    to the database although it was already saved in it.

.. class:: sunpy.database.EntryAlreadyStarredError(database_entry)

    This exception is raised if a database entry is marked as starred
    using :meth:`Database.star` although it was already starred before
    this operation.

.. class:: sunpy.database.EntryAlreadyUnstarredError(database_entry)

    This exception is raised if the star mark from a database entry is
    attempted to be removed although the entry is not starred.

.. class:: sunpy.database.EntryNotFoundError(entry_id)

    This exception is raised if a database entry cannot be found by its
    unique ID.

.. class:: sunpy.database.TagAlreadyAssignedError(database_entry, tag_name)

    This exception is raised if it is attempted to assign a tag to a
    database entry but the database entry already has this tag assigned.

.. class:: sunpy.database.NoSuchTagError(tag_name)

    This exception is raised if a tag cannot be found in a database by its
    name.

Submodules
----------

.. automodapi:: sunpy.database.tables
    :headings: ^+
    :no-inheritance-diagram:

.. automodapi:: sunpy.database.caching
    :headings: ^+
    :no-inheritance-diagram:

.. automodapi:: sunpy.database.commands
    :headings: ^+
    :no-inheritance-diagram:

.. automodapi:: sunpy.database.attrs
    :headings: ^+
    :no-inheritance-diagram:

