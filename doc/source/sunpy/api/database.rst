SunPy Database
==============
..
    TODO:
        - WHY was this package developed?
        - WHAT are the use-cases?
        - WHICH methods can be undone?
        - WHICH methods are the most important ones to know?

.. automodapi:: sunpy.database

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

sunpy.database.tables Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: sunpy.database.tables.DatabaseEntry
   :members:

.. autoclass:: sunpy.database.tables.FitsHeaderEntry
   :members:

.. autoclass:: sunpy.database.tables.FitsKeyComment
   :members:

.. autoclass:: sunpy.database.tables.Tag
   :members:

.. autoclass:: sunpy.database.tables.JSONDump
   :members:

utility functions
+++++++++++++++++
.. autofunction:: sunpy.database.tables.entries_from_query_result

.. autofunction:: sunpy.database.tables.entries_from_file

.. autofunction:: sunpy.database.tables.entries_from_dir

.. autofunction:: sunpy.database.tables.display_entries

Exceptions
++++++++++
.. class:: sunpy.database.tables.WaveunitNotFoundError(obj)

    This exception is raised if a wavelength unit cannot be found in a FITS
    header or in a VSO query result block.

.. class:: sunpy.database.tables.WaveunitNotConvertibleError(waveunit)

    This exception is raised if a wavelength cannot be converted to an
    astropy.units.Unit instance.

.. automodapi:: sunpy.database.caching
    :headings: "^+"

sunpy.database.commands Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: sunpy.database.commands
   :exclude-members:
        EmptyCommandStackError, NoSuchEntryError, NonRemovableTagError
   :members:
   :show-inheritance:

Exceptions
++++++++++
.. class:: sunpy.database.commands.EmptyCommandStackError

    This exception is raised if it is attempted to pop from a command
    stack even though it is empty.

.. class:: sunpy.database.commands.NoSuchEntryError(database_entry)

    This exception is raised if it is attempted to remove an entry even
    though it does not exist in the database.

.. class:: sunpy.database.commands.NonRemovableTagError(database_entry, tag)

    This exception is raised if it is attempted to remove a tag from a
    database entry even though it is not saved in this entry.

sunpy.database.attrs Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. module:: sunpy.database.attrs

.. class:: sunpy.database.attrs.Starred

.. class:: sunpy.database.attrs.Tag(tagname)

.. class:: sunpy.database.attrs.Path(value)

.. class:: sunpy.database.attrs.DownloadTime(start, end)

.. class:: sunpy.database.attrs.FitsHeaderEntry(key, value)
