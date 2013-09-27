TODO
====

Documentation
-------------

Reference
~~~~~~~~~
- Database: explain the purpose and meaning of this class

- DatabaseEntry:

  - explain the purpose and meaning of this class

  - say that instances of this class are usually not created manually but
    by other methods such as Database.add_from_path or
    Database.add_from_vso_query_result

- document the tables Tag, FitsHeaderEntry, FitsKeyComment, JSONDump

- attrs module: document the following classes:

  - Starred

  - Tag

  - Path

  - DownloadTime

  - FitsHeaderEntry

- attrs: document that for vso.attrs.Time, the `near` parameter is
  ignored!

- commands.CommandManager: update docstrings, especially paying attention
  to using lists instead of commands as parameters!

- document property cache_size, cache_maxsize

- document Database.__getitem__

- improve doc of commands.DatabaseOperation (parameters of __init__)

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

Important
---------
- support serialize.load_{query,entry} with datetime format!

- serialize: support indent parameter in dump_query

- test_database: remove ``database_without_tables``

- Database.{add_from_file,add_from_path},
  entries_from_file, entries_from_dir: support *ignoring* entries without
  waveunit

- display_entries: support and test output of FITS header entries

- adopt parameter names in Database methods to last changes in the
  documentation!!!

- support (un-)pickling of instances from the tables module

Low priority
---------
- support saving entries by HEK query result

- check if catching the exception InvalidRequestError is really sufficient
  in AddEntry.__call__, AddEntry.undo, RemoveEntry.__call__

- support using the progressbar for adding many entries to the database
  → Or do not support it directly. Rather, add a section in the docs "Tips
  & Tricks" or something and show a snippet there

- support TAR and TAR-GZ

Unclear / Undecided
-------------------
- Should the FITS Header be normalized before being saved (the keys are
  usually in UPPERCASE, makes it worse to search for)?

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
