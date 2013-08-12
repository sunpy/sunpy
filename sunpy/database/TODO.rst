TODO
====

Documentation
-------------

Reference
~~~~~~~~~
- commands: use :exc:`...` to refer to exceptions

- Database: explain the purpose and meaning of this class

- DatabaseEntry:

  - explain the purpose and meaning of this class

  - say that instances of this class are usually not created manually but
    by other methods such as Database.add_from_path or
    Database.add_from_vso_query_result

- attrs module: document the following classes:
  
  - Starred
    
  - Tag
    
  - Path
    
  - DownloadTime
    
  - FitsHeaderEntry

- document in more methods and functions which exceptions may be raised

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
- test undoing and redoing the following methods of Database:

  - tag

  - star

  - unstar

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
- make BaseCache.dict private → BaseCache._dict

- do not save the waveunit. rather, use only nm -> only possible if
  PR #522 is merged

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
- write a contextmanager to disable undo/redo functionality for a block of
  operations -> commands module

- support (un-)pickling of instances from the tables module

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

HEK Notes
---------
important HEK result keys:

    - u'obs_wavelunit'

    - u'obs_instrument'

    - u'event_starttime'

    - u'event_endtime'

    - u'obs_observatory'?

important attributes: vso_time, vso_instrument

relevant function: translate_results_to_query from the hek2vso package
