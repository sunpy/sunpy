TODO
====

Documentation
-------------

Reference
~~~~~~~~~
- Database: explain the purpose and meaning of this class

- DatabaseEntry: explain the purpose and meaning of this class

- document properties cache_size, cache_maxsize

- document Database.__getitem__

Testing
-------
- what to do if a database operation is attempted but the table hasn't
  been created yet? Throw a custom exception or create the required
  table(s) silently? Currently, the exception OperationalError from
  sqlalchemy is raised in such a case.

Important
---------
- support serialize.load_{query,entry} with datetime format!

- serialize: support indent parameter in dump_query

- Database.{add_from_file,add_from_path},
  entries_from_file, entries_from_dir: support *ignoring* entries without
  waveunit

- display_entries: support and test output of FITS header entries

Low priority
---------
- support saving entries by HEK query result

- support TAR and TAR-GZ
