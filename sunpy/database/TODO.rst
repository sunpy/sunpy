TODO
====

Documentation
-------------
- document property cache_size, cache_maxsize

- document commands.EmptyCommandStackError, commands.NoSuchEntryError

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
- trivial: new Database.clear method to remove all entries in one run
  (again, undo/redo should be intuitevely possible)

- do not save the waveunit. rather, use only nm

- support the VSO attribute Time for querying the database

- support multiple HDUs per database entry

- support saving entries by HEK query result

- __repr__ for the classes in commands module

- display_entries: support and test output of FITS header entries

- use new-style string formatting everywhere

- adopt parameter names in Database methods to last changes in the
  documentation!!!

- check if catching the exception InvalidRequestError is really sufficient
  in AddEntry.__call__, AddEntry.undo, RemoveEntry.__call__

Low priority
---------

- Database.{__contains__, __iter__, __len__}: return values of the cache
  instead of sending a query to the database
  PRO: better speed performance
  CONTRA: needs explicit "commit" calls after each database operation or
  before each __contains__, __iter__, __len__ call

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
