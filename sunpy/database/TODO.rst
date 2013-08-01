TODO
====

Documentation
-------------

- document the package database itself!!! → show how to use the Database
  class and the undo and redo methods

- document Database.__getitem__

- document the caching package

- document caching.LRUCache and caching.LFUCache

- improve doc of commands.DatabaseOperation (parameters of __init__)

- do not autodocument the inherited methods of caching.BaseCache

- document the specific commands AddEntry, EditEntry, RemoveEntry

Testing
-------
- test undoing and redoing the latest methods in Database!!! (e.g.
  tagging, removing tags, ...)

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

Important
---------
- do not save the waveunit. rather, use only nm

- display_entries: support and test output of FITS header entries

- support the VSO attribute Time for querying the database

- adopt parameter names in Database methods to last changes in the
  documentation!!!

- which table fields are unique?

- check if catching the exception InvalidRequestError is really sufficient
  in AddEntry.__call__, AddEntry.undo, RemoveEntry.__call__

Low priority
---------

- support using the progressbar for adding many entries to the database
  → Or do not support it directly. Rather, add a section in the docs "Tips
  & Tricks" or something and show a snippet there

- support TAR and TAR-GZ

Unclear / Undecided
-------------------

- Should the FITS Header be normalized before being saved (the keys are
  usually in UPPERCASE, makes it worse to search for)?

Interfacing with VSO
--------------------
- general problem: db query cannot tell whether the same to the VSO might
  fetch more data or the same data

  → solution: remember queries, compare
