TODO
====

Documentation
-------------

- document the package database itself!!! → show how to use the Database
  class and the undo and redo methods

- document entries_from_query_result

- document the caching package

- document caching.LRUCache and caching.LFUCache

- improve doc of commands.DatabaseOperation (parameters of __init__)

- do not autodocument the inherited methods of caching.BaseCache

- document the specific commands AddEntry, EditEntry, RemoveEntry

Testing
-------

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

- use a main DatabaseException class for all custom exceptions that are
  raised in this package → also a main sunpy.SunPyException

- support changing the cache size from Database → Database.cache_size

Important
---------

- adopt parameter names in Database methods to last changes in the
  documentation!!!

- support removing tags → see sqlalchemy doc to see what to pay attention
  to when removing entries in a many-to-many relationship

- which table fields are unique?

- check if catching the exception InvalidRequestError is really sufficient
  in AddEntry.__call__, AddEntry.undo, RemoveEntry.__call__

Long-term
---------

- support querying via attributes (see vso)

Unclear / Undecided
-------------------

- Should the FITS Header be normalized before being saved (the keys are
  usually in UPPERCASE, makes it worse to search for)?

- how should the database-vso connection be implemented exactly?
