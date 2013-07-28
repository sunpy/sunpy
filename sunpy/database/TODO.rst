TODO
====

Documentation
-------------

- document the package database itself!!! → show how to use the Database
  class and the undo and redo methods

- document the caching package

- document caching.LRUCache and caching.LFUCache

- improve doc of commands.DatabaseOperation (parameters of __init__)

- do not autodocument the inherited methods of caching.BaseCache

- document the specific commands AddEntry, EditEntry, RemoveEntry

Testing
-------
- test adding entries that are already saved in the database → should
  throw an exception

- what to do if a database operation is attempted but the table hasn't
  been created yet? Throw a custom exception or create the required
  table(s) silently? Currently, the exception OperationalError from
  sqlalchemy is raised in such a case.

- test commands.CommandManager
  → test undo/redo in combination with the caching strategy
    (integration test, wohoo)

- test starring and tagging entries that haven't been added to the
  database!

Simple
------
- support changing the cache size from Database → Database.cache_size

Important
---------
- Database.query / attrs module:
  
  - support VSO attributes

  - support more attributes:

      - Path (str, compare Tag attribute)

      - DownloadTime (time range, see vso.attrs.Time)

      - FitsHeaderEntry (custom attribute storing key and value)

- remove Database.get_entry_by_id

- adopt parameter names in Database methods to last changes in the
  documentation!!!

- support removing tags → see sqlalchemy doc to see what to pay attention
  to when removing entries in a many-to-many relationship

- which table fields are unique?

- check if catching the exception InvalidRequestError is really sufficient
  in AddEntry.__call__, AddEntry.undo, RemoveEntry.__call__

Long-term
---------

- support using the progressbar for adding many entries to the database
  → Or do not support it directly. Rather, add a section in the docs "Tips
  & Tricks" or something and show a snippet there

- support TAR and TAR-GZ

Unclear / Undecided
-------------------

- Should the FITS Header be normalized before being saved (the keys are
  usually in UPPERCASE, makes it worse to search for)?
