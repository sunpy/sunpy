.. _database_guide:

--------------------------
Using the database package
--------------------------
.. currentmodule:: sunpy.database

The database package offers the possibility to save retrieved data (e.g.
via the :mod:`sunpy.net` package) onto a local or remote database. The
database may be a single file located on a local hard drive (if a SQLite
database is used) or a local or remote database server (see the SQLAlchemy
documentation for a `list of supported databases
<http://docs.sqlalchemy.org/en/stable/dialects/>`_)
This makes it possible to fetch required data from the local database
instead of downloading it again from a remote server.

The package :mod:`sunpy.database` was developed as part of
`Google Summer of Code (GSOC) 2013
<https://www.google-melange.com/archive/gsoc/2013>`_.

1. Connecting and initializing the database
-------------------------------------------
To start a connection to an existing or a new database, create
a :class:`Database` object:

    >>> from sunpy.database import Database
    >>> database = Database('sqlite:///sunpydata.sqlite')

The database object in our example above connects to a new SQLite database with
the file name "sunpydata.sqlite" in the current directory.

The first parameter of :class:`Database` receives one mandatory argument:
a URL which describes how to connect to the database. The supported
format of this URL is described by the documentation of
:func:`sqlalchemy.create_engine`.

Note that a connection is only established when it's really needed,
i.e. if some query is sent to the database to read from it. Transactions
can also be committed explicitly using the :meth:`Database.commit` method.

.. warning::

    If you are using :class:`Database` objects in an interactive Python
    session you must not forget to call the :meth:`Database.commit` method
    on them explicitly before quitting the Python session! Otherwise, all
    changes on the altered databases are lost!

.. note::

    You can set the default database url in the sunpy config file under the
    'database' section. See :ref:`customizing-sunpy` for information on the
    config file. A database section might look like this::

        [database]
        url = sqlite:////home/user/sunpy/my_database.sqlite


2. Adding new entries
---------------------
Each entry in a database is an instance of the class
:class:`tables.DatabaseEntry` with the following attributes:

====================== ===================================================
      Attribute                            Description
====================== ===================================================
id                     A unique ID number. By default it is None, but
                       automatically set to the maximum number plus one
                       when an entry is added to the database.
source                 The source is the name of an observatory or the
                       name of a network of observatories.
provider               The name of the server which provides the retrieved
                       data.
physobs                A physical observable identifier used by VSO.
fileid                 The file ID is a string defined by the data
                       provider that should point to a specific data
                       product. The association of fileid to the specific
                       data may change sometime, if the fileid always
                       points to the latest calibrated data.
observation_time_start The date and time when the observation of the data
                       started.
observation_time_end   The date and time when the observation of the data
                       ended.
instrument             The instrument which was used to observe the data.
size                   The size of the data in kilobytes (-1 if unknown).
wavemin                The value of the measured wave length.
wavemax                This is the same value as ``wavemin``. The value is
                       stored twice, because each
                       ``suds.sudsobject.QueryResponseBlock`` which is
                       used by the vso package contains both these values.
path                   A local file path where the according FITS file is
                       saved.
download_time          The date and time when the files connected to a
                       query have been downloaded. Note: this is not the
                       date and time when this entry has been added to a
                       database!
starred                Entries can be starred to mark them. By default,
                       this value is False.
fits_header_entries    A list of :class:`tables.FitsHeaderEntry` instances.
fits_key_comments      A list of :class:`tables.FitsKeyComment` instances.
tags                   A list of :class:`tables.Tag` instances.
====================== ===================================================

* The ``id`` attribute is automatically set if an entry is added to a database.

* The attributes ``source``, ``provider``, ``physobs``, ``fileid``,
  ``observation_time_start``, ``observation_time_end``, ``instrument``,  ``size``,
  ``wavemin``, and ``wavemax`` are set by methods which use the VSO interface. In
  particular, these are :meth:`Database.add_from_vso_query_result`,
  :meth:`Database.download` and   possibly :meth:`Database.fetch`.

* The attributes ``path`` and ``download_time`` are set by the method
  :meth:`Database.download` and also possibly by :meth:`Database.fetch`.

* ``starred`` is set or changed via :meth:`Database.star` or :meth:`unstar`.

* ``tags`` is set via :meth:`Database.tag` or :meth:`Database.remove_tag`.

* The attribute ``fits_header_entries`` is set by the methods
  :meth:`Database.download`, :meth:`Database.add_from_dir`, and
  :meth:`Database.add_from_file`.

2.1 Adding entries from one FITS file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The method :meth:`Database.add_from_file` receives one positional argument
(either a string or a file-like object) which is used to add at least one
new entry from the given FITS file to the database. Why "at least one" and
not "exactly one"? The reason is that each database entry does not
represent one file, but one FITS header from a file. That means if you pass
a file which has 5 FITS headers in it, 5 entries will be added to the
database. The file in the following example
(``sunpy.data.sample.AIA_171_IMAGE``) has only one FITS header, so just one
entry is added to the database.

.. note::

  If you are working with Hinode/SOT files you may notice that for
  each file you get two entries, one which refers to the observation and
  another that contains some (useless - as `discussed
  <http://listmgr.cv.nrao.edu/pipermail/fitsbits/2007-August/001923.html>`_
  in the `fitsbits mailing list
  <http://listmgr.cv.nrao.edu/mailman/listinfo/fitsbits>`_)
  telemetry data.

:meth:`Database.add_from_file` saves the value of ``path`` by either simply
passing on the value of the received argument (if it was a string)
or by reading the value of ``file.name`` where ``file`` is the passed argument.
If the path cannot be determined, it stays as ``None`` (the default value).

``wavemin`` and `wavemax`` are only set if the wavelength unit
of the FITS file can be found out or if the ``default_waveunit`` attribute of
the database object is set. These values are then
used to convert from the used units to nanometers. The rationale behind
this behaviour is that it makes querying for wavelengths more flexible.

``instrument`` is set by looking up the
FITS header key *INSTRUME*. ``observation_time_start`` is set
by searching for the FITS header key *DATE-OBS* or *DATE_OBS*.
Analogously, ``observation_time_end`` is set by searching for *DATE-END* or
*DATE_END*. Finally, the whole FITS header is stored in the attribute
``fits_header_entries`` as a list of :class:`tables.FitsHeaderEntry`
instances, and FITS comments are stored in the attribute
``fits_key_comments`` which is a list of :class:`tables.FitsKeyComment`
instances.

Using the function `len` on a :class:`Database` object returns the
number of saved database entries. To get the first entry of the database,
``database[0]`` is used (the ID number of the entries does not matter,
``database[0]`` always returns the oldest saved entry of the database). If
the database is empty, this expression raises an :exc:`IndexError`.
In section 3, more advanced formats of the slicing syntax are introduced.

    >>> import sunpy.data.sample
    >>> database.add_from_file(sunpy.data.sample.AIA_171_IMAGE)
    >>> len(database)
    1
    >>> entry = database[0]
    >>> entry.path == sunpy.data.sample.AIA_171_IMAGE
    True
    >>> entry.wavemin, entry.wavemax
    (17.1, 17.1)
    >>> entry.instrument
    'AIA_3'
    >>> entry.observation_time_start, entry.observation_time_end
    (datetime.datetime(2011, 6, 7, 6, 33, 2, 770000), None)
    >>> len(entry.fits_header_entries)
    191
    >>> for fits_header_entry in entry.fits_header_entries[:10]:
    ...     print('{entry.key}\n\t{entry.value}'.format(entry=fits_header_entry))   # doctest: +NORMALIZE_WHITESPACE
    SIMPLE
    	True
    BITPIX
    	-32
    NAXIS
    	2
    NAXIS1
    	1024
    NAXIS2
    	1024
    BLD_VERS
    	V5R12X
    LVL_NUM
    	1.5
    T_REC
    	2011-06-07T06:33:03Z
    TRECSTEP
    	1.0
    TRECEPOC
    	1977.01.01_00:00:00_TAI

    >>> for fits_key_comment in entry.fits_key_comments:
    ...     print('{comment.key}\n\t{comment.value}'.format(comment=fits_key_comment))   # doctest: +NORMALIZE_WHITESPACE
    SIMPLE
    	conforms to FITS standard
    BITPIX
    	array data type
    NAXIS
    	number of array dimensions


2.2 Adding entries from a directory of FITS files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Adding all FITS files from a directory works by calling the method
:meth:`Database.add_from_dir` and passing the desired directory to it. By
setting the keyword argument ``ignore_already_added`` to ``True``, no
exception is raised if it is attempted to add an already existing entry

In the following example case, setting this parameter is required because the
file ``sunpy.data.sample.AIA_171_IMAGE`` was already added, which is located
in the directory ``sampledata_dir``

    >>> from sunpy import config
    >>> sampledata_dir = config.get("downloads", "sample_dir")
    >>> database.default_waveunit = 'angstrom'
    >>> database.add_from_dir(sampledata_dir, ignore_already_added=True,
    ...                       time_string_parse_format="%d/%m/%Y")
    >>> len(database)
    58

2.3 Adding entries using the VSO interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

2.3.1 Adding entries from a VSO query result
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The number of database entries that will be added from a VSO query result
is equal to the value of ``len(qr)`` in the following code snippet. Note that
:meth:`Database.add_from_vso_query_result` does not download any files,
though. If you want to add new entries using the VSO and also want to
download files at the same time, take a look at the following two
sections.

    >>> from sunpy.net import vso
    >>> client = vso.VSOClient()  # doctest: +REMOTE_DATA
    >>> qr = client.search(
    ...     vso.attrs.Time('2011-05-08', '2011-05-08 00:00:05'),
    ...     vso.attrs.Instrument('AIA'))  # doctest: +REMOTE_DATA
    >>> len(qr)  # doctest: +REMOTE_DATA
    4
    >>> database.add_from_vso_query_result(qr)  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    62

2.3.2 "Clever" Fetching
^^^^^^^^^^^^^^^^^^^^^^^

The method :meth:`Database.fetch` checks if the given query has already been
used once to add entries to the database. Otherwise, the query is used to
download and add new data. The :meth:`Database.fetch` method also accepts an
optional keyword argument ``path`` which is passed as-is to
:meth:`sunpy.net.vso.VSOClient.get` and determines the value of the ``path``
attribute of each entry.

Note that the number of entries that will be added depends on the total number
of FITS headers, **not** the number of records in the query.

In the next code snippet new data is downloaded as this query
has not been downloaded before.

    >>> entries = database.fetch(
    ...     vso.attrs.Time('2012-08-05', '2012-08-05 00:00:05'),
    ...     vso.attrs.Instrument('AIA'))  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    66

Any more identical fetch calls don't add any new entries to the database
because they have already been downloaded.

    >>> entries = database.fetch(
    ...     vso.attrs.Time('2012-08-05', '2012-08-05 00:00:05'),
    ...     vso.attrs.Instrument('AIA'))  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    66

However, this different fetch call downloads new files because there is a
new date range whose files have not been downloaded yet.

    >>> entries = database.fetch(
    ...     vso.attrs.Time('2013-08-05', '2013-08-05 00:00:05'),
    ...     vso.attrs.Instrument('AIA'))  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    70

The caching also ensures that when queries have some results in
common, files for the common results will not be downloaded again. In the
following example, the query is new, but all of it's files have already been
downloaded. This means no new files are downloaded.

    >>> entries = database.fetch(
    ...     vso.attrs.Time('2012-08-05 00:00:00', '2012-08-05 00:00:01'),
    ...     vso.attrs.Instrument('AIA'))  # doctest: +REMOTE_DATA
    >>> len(database)  # doctest: +REMOTE_DATA
    70

2.4 Adding entries manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Although usually not required, it is also possible to add database entries
by specifying the parameters manually. To do so pass the
values as keyword arguments to :class:`tables.DatabaseEntry` as follows:

    >>> from sunpy.database.tables import DatabaseEntry
    >>> entry = DatabaseEntry(instrument='EIT', wavemin=25.0)
    >>> database.add(entry)
    >>> entry in database
    False
    >>> database.commit()
    >>> entry in database
    True
    >>> len(database) # doctest: +REMOTE_DATA
    71

Note that the `in` operator only works as expected after the
:meth:`Database.commit` method has been called!

3. Displaying entries in a table
--------------------------------
In the previous code snippets 71 entries have been added,
all of them saving a lot of data. To display the database in a table format
there is a helper function. :func:`tables.display_entries` takes two
arguments: the first one is an iterator of :class:`tables.DatabaseEntry`
instances. Remember that an instance of :class:`Database` yields instances
of :class:`tables.DatabaseEntry` instances, so you can simply pass a
database object. The second argument is an iterable of the resulting
columns in the table to be displayed. Each string in this iterable is used
to access the entry's attribute values. In the following example, the
values of ``entry.id``, ``entry.observation_time_start``,
``entry.observation_time_end``, ``entry.instrument``, ``entry.wavemin``,
and ``entry.wavemax`` are displayed (where ``entry`` stands for the
respective database entry). Note that *N/A* is displayed if the value
cannot be found or is not set.

    >>> from sunpy.database.tables import display_entries
    >>> print(display_entries(database,
    ...                       ['id', 'observation_time_start', 'observation_time_end',
    ...                        'instrument', 'wavemin', 'wavemax']))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start ...      wavemin            wavemax
    --- ---------------------- ... ------------------ ------------------
      1    2011-06-07 06:33:02 ...               17.1               17.1
      2    2011-06-07 06:33:01 ... 13.100000000000001 13.100000000000001
      3    2011-06-07 06:33:02 ...               17.1               17.1
      4    2011-06-07 06:33:02 ...               21.1               21.1
      5    2011-06-07 06:33:03 ...               33.5               33.5
      6    2011-06-07 06:33:05 ...                9.4                9.4
      7    2011-06-07 06:33:05 ...              160.0              160.0
      8    2011-06-07 06:33:07 ...               19.3               19.3
      9    2011-06-07 06:33:07 ...               19.3               19.3
     10    2011-06-07 06:39:31 ...               19.3               19.3
    ...                    ... ...                ...                ...
     61    2011-05-08 00:00:02 ...                9.4                9.4
     62    2011-05-08 00:00:03 ...               33.5               33.5
     63    2012-08-05 00:00:01 ...                9.4                9.4
     64    2012-08-05 00:00:01 ...                9.4                9.4
     65    2012-08-05 00:00:02 ...               33.5               33.5
     66    2012-08-05 00:00:02 ...               33.5               33.5
     67    2013-08-05 00:00:01 ...                9.4                9.4
     68    2013-08-05 00:00:01 ...                9.4                9.4
     69    2013-08-05 00:00:02 ...               33.5               33.5
     70    2013-08-05 00:00:02 ...               33.5               33.5
     71                    N/A ...               25.0                N/A
    Length = 71 rows

The index operator can be used to access certain single database entries like
``database[5]``  to get the 5th entry or  ``database[-2]``
to get the second-latest entry. It is also possible to
use more advanced slicing syntax;
see https://docs.scipy.org/doc/numpy/reference/arrays.indexing.html
for more information.

    >>> print(display_entries(database[9::10],
    ...                       ['id', 'observation_time_start', 'observation_time_end',
    ...                        'instrument', 'wavemin', 'wavemax']))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start observation_time_end instrument wavemin wavemax
    --- ---------------------- -------------------- ---------- ------- -------
     10    2011-06-07 06:39:31                  N/A      AIA_2    19.3    19.3
     20                    N/A                  N/A        N/A     N/A     N/A
     30    2011-06-06 23:59:55  2011-06-08 00:00:05        GBM     N/A     N/A
     40                    N/A                  N/A        N/A     N/A     N/A
     50                    N/A                  N/A        N/A     N/A     N/A
     60    2011-05-08 00:00:00  2011-05-08 00:00:01        AIA    21.1    21.1
     70    2013-08-05 00:00:02  2013-08-05 00:00:03        AIA    33.5    33.5

4. Removing entries
-------------------
`database.remove()` can be used to remove database entries from the SunPy
database. It takes a ``tables.DatabaseEntry`` object as argument.

For example, imagine we want to only have database entries which have an
observation time saved. To remove all the entries where the value of both
``observation_time_start`` and ``observation_time_end`` is ``None``,
simply iterate over the database and the :meth:`Database.remove`
method to remove those where there is no time set:

    >>> for database_entry in database:
    ...     if database_entry.observation_time_start is None and database_entry.observation_time_end is None:
    ...         database.remove(database_entry)
    ...
    >>> len(database) # doctest: +REMOTE_DATA
    38
    >>> print(display_entries(database,
    ...                       ['id', 'observation_time_start', 'observation_time_end',
    ...                        'instrument', 'wavemin', 'wavemax']))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start ...      wavemin            wavemax
    --- ---------------------- ... ------------------ ------------------
      1    2011-06-07 06:33:02 ...               17.1               17.1
      2    2011-06-07 06:33:01 ... 13.100000000000001 13.100000000000001
      3    2011-06-07 06:33:02 ...               17.1               17.1
      4    2011-06-07 06:33:02 ...               21.1               21.1
      5    2011-06-07 06:33:03 ...               33.5               33.5
      6    2011-06-07 06:33:05 ...                9.4                9.4
      7    2011-06-07 06:33:05 ...              160.0              160.0
      8    2011-06-07 06:33:07 ...               19.3               19.3
      9    2011-06-07 06:33:07 ...               19.3               19.3
     10    2011-06-07 06:39:31 ...               19.3               19.3
    ...                    ... ...                ...                ...
     60    2011-05-08 00:00:00 ...               21.1               21.1
     61    2011-05-08 00:00:02 ...                9.4                9.4
     62    2011-05-08 00:00:03 ...               33.5               33.5
     63    2012-08-05 00:00:01 ...                9.4                9.4
     64    2012-08-05 00:00:01 ...                9.4                9.4
     65    2012-08-05 00:00:02 ...               33.5               33.5
     66    2012-08-05 00:00:02 ...               33.5               33.5
     67    2013-08-05 00:00:01 ...                9.4                9.4
     68    2013-08-05 00:00:01 ...                9.4                9.4
     69    2013-08-05 00:00:02 ...               33.5               33.5
     70    2013-08-05 00:00:02 ...               33.5               33.5
    Length = 38 rows


5. Editing entries
------------------

5.1 Starring and unstarring entries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The database package supports marking certain entries as "starred" using the
:meth:`Database.star` method. For example, to star
all values that have a wavelength of 20nm or higher:

    >>> for database_entry in database:
    ...     if database_entry.wavemin and database_entry.wavemin > 20:
    ...         database.star(database_entry)
    >>> print(display_entries(
    ...     filter(lambda entry: entry.starred, database),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax']))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start observation_time_end instrument wavemin wavemax
    --- ---------------------- -------------------- ---------- ------- -------
      4    2011-06-07 06:33:02                  N/A      AIA_2    21.1    21.1
      5    2011-06-07 06:33:03                  N/A      AIA_1    33.5    33.5
      7    2011-06-07 06:33:05                  N/A      AIA_3   160.0   160.0
     60    2011-05-08 00:00:00  2011-05-08 00:00:01        AIA    21.1    21.1
     62    2011-05-08 00:00:03  2011-05-08 00:00:04        AIA    33.5    33.5
     65    2012-08-05 00:00:02  2012-08-05 00:00:03        AIA    33.5    33.5
     66    2012-08-05 00:00:02  2012-08-05 00:00:03        AIA    33.5    33.5
     69    2013-08-05 00:00:02  2013-08-05 00:00:03        AIA    33.5    33.5
     70    2013-08-05 00:00:02  2013-08-05 00:00:03        AIA    33.5    33.5

To remove the star from these entries, the :meth:`Database.unstar` method
works the same way.

5.2 Setting and removing tags
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To add some more information by assigning tags to entries, the database package
also supports *tags*. The ``tags`` property of a database object holds all tags
that are saved in the database. For example, to assign the tag *spring* to all
entries that have been observed in March, April, or May on any day in any
year:

    >>> for database_entry in database:
    ...     if database_entry.observation_time_start.month in [3,4,5]:
    ...         database.tag(database_entry, 'spring')
    >>> database.tags
    [<Tag(name 'spring')>]
    >>> spring = database.tags[0]
    >>> print(display_entries(
    ...     filter(lambda entry: spring in entry.tags, database),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax']))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start observation_time_end instrument wavemin wavemax
    --- ---------------------- -------------------- ---------- ------- -------
     22    2014-04-09 06:00:12                  N/A      AIA_3    17.1    17.1
     59    2011-05-08 00:00:00  2011-05-08 00:00:01        AIA    17.1    17.1
     60    2011-05-08 00:00:00  2011-05-08 00:00:01        AIA    21.1    21.1
     61    2011-05-08 00:00:02  2011-05-08 00:00:03        AIA     9.4     9.4
     62    2011-05-08 00:00:03  2011-05-08 00:00:04        AIA    33.5    33.5

5.3 Manually changing attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Attributes for entries in the database can be manually edited. The
:meth:`Database.edit` method receives the database entry to be edited and any
number of keyword arguments to describe which values to change and how. For
example, to change the time entires that are ``None`` to the start of the
observation time (because it's possible, not because it is accurate!):

    >>> for database_entry in database:
    ...     if database_entry.observation_time_end is None:
    ...         database.edit(database_entry, observation_time_end=database_entry.observation_time_start)
    ...
    >>> print(display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax']))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start ...      wavemin            wavemax
    --- ---------------------- ... ------------------ ------------------
      1    2011-06-07 06:33:02 ...               17.1               17.1
      2    2011-06-07 06:33:01 ... 13.100000000000001 13.100000000000001
      3    2011-06-07 06:33:02 ...               17.1               17.1
      4    2011-06-07 06:33:02 ...               21.1               21.1
      5    2011-06-07 06:33:03 ...               33.5               33.5
      6    2011-06-07 06:33:05 ...                9.4                9.4
      7    2011-06-07 06:33:05 ...              160.0              160.0
      8    2011-06-07 06:33:07 ...               19.3               19.3
      9    2011-06-07 06:33:07 ...               19.3               19.3
     10    2011-06-07 06:39:31 ...               19.3               19.3
    ...                    ... ...                ...                ...
     60    2011-05-08 00:00:00 ...               21.1               21.1
     61    2011-05-08 00:00:02 ...                9.4                9.4
     62    2011-05-08 00:00:03 ...               33.5               33.5
     63    2012-08-05 00:00:01 ...                9.4                9.4
     64    2012-08-05 00:00:01 ...                9.4                9.4
     65    2012-08-05 00:00:02 ...               33.5               33.5
     66    2012-08-05 00:00:02 ...               33.5               33.5
     67    2013-08-05 00:00:01 ...                9.4                9.4
     68    2013-08-05 00:00:01 ...                9.4                9.4
     69    2013-08-05 00:00:02 ...               33.5               33.5
     70    2013-08-05 00:00:02 ...               33.5               33.5
    Length = 38 rows

It is possible to use
``database_entry.observation_time_end = database_entry.observation_time_start``",
but it has one major disadvantage: this operation cannot be undone.
See section 6 to see how undoing and redoing works.

6. Undoing and redoing operations
---------------------------------
A very handy feature of the database package is that every operation that
changes the database in some way can be reverted. The Database class has
the methods :meth:`Database.undo` and :meth:`Database.redo` to undo
and redo the last n commands, respectively.

.. warning::

    The undo and redo history are only saved in-memory! That means in
    particular, that if you work on a Database object in an interactive
    Python session and quit this session, the undo and redo history are
    lost.

In the following snippet, the operations from the sections 5.3, 5.2, and
5.1 are undone. Note the following changes: there are no longer any tags
anymore saved in the database, there is no entry which is starred and
the entries with no end of observation time are back.

    >>> database.undo(4)  # undo the edits from 5.3 (4 records have been affected)
    >>> print(display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax', 'tags', 'starred']))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start observation_time_end ...  tags  starred
    --- ---------------------- -------------------- ... ------ -------
      1    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A      No
      2    2011-06-07 06:33:01  2011-06-07 06:33:01 ...    N/A      No
      3    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A      No
      4    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A     Yes
      5    2011-06-07 06:33:03  2011-06-07 06:33:03 ...    N/A     Yes
      6    2011-06-07 06:33:05  2011-06-07 06:33:05 ...    N/A      No
      7    2011-06-07 06:33:05  2011-06-07 06:33:05 ...    N/A     Yes
      8    2011-06-07 06:33:07  2011-06-07 06:33:07 ...    N/A      No
      9    2011-06-07 06:33:07  2011-06-07 06:33:07 ...    N/A      No
     10    2011-06-07 06:39:31  2011-06-07 06:39:31 ...    N/A      No
    ...                    ...                  ... ...    ...     ...
     60    2011-05-08 00:00:00  2011-05-08 00:00:01 ... spring     Yes
     61    2011-05-08 00:00:02  2011-05-08 00:00:03 ... spring      No
     62    2011-05-08 00:00:03  2011-05-08 00:00:04 ... spring     Yes
     63    2012-08-05 00:00:01  2012-08-05 00:00:02 ...    N/A      No
     64    2012-08-05 00:00:01  2012-08-05 00:00:02 ...    N/A      No
     65    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    N/A     Yes
     66    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    N/A     Yes
     67    2013-08-05 00:00:01  2013-08-05 00:00:02 ...    N/A      No
     68    2013-08-05 00:00:01  2013-08-05 00:00:02 ...    N/A      No
     69    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    N/A     Yes
     70    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    N/A     Yes
    Length = 38 rows

The :meth:`~Database.redo` method reverts the last n operations that have been
undone. If not that many operations can be redone (i.e. any number greater than
14 in this example), an exception is raised. The redo call
reverts the original state: the tags have appeared again and all entries with a
wavelength >20nm are starred again, and there are no entries with no
stored end of observation time.

    >>> database.redo(1)  # redo all undone operations
    >>> print(display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax', 'tags', 'starred']))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start observation_time_end ...  tags  starred
    --- ---------------------- -------------------- ... ------ -------
      1    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A      No
      2    2011-06-07 06:33:01  2011-06-07 06:33:01 ...    N/A      No
      3    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A      No
      4    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    N/A     Yes
      5    2011-06-07 06:33:03  2011-06-07 06:33:03 ...    N/A     Yes
      6    2011-06-07 06:33:05  2011-06-07 06:33:05 ...    N/A      No
      7    2011-06-07 06:33:05  2011-06-07 06:33:05 ...    N/A     Yes
      8    2011-06-07 06:33:07  2011-06-07 06:33:07 ...    N/A      No
      9    2011-06-07 06:33:07  2011-06-07 06:33:07 ...    N/A      No
     10    2011-06-07 06:39:31  2011-06-07 06:39:31 ...    N/A      No
    ...                    ...                  ... ...    ...     ...
     60    2011-05-08 00:00:00  2011-05-08 00:00:01 ... spring     Yes
     61    2011-05-08 00:00:02  2011-05-08 00:00:03 ... spring      No
     62    2011-05-08 00:00:03  2011-05-08 00:00:04 ... spring     Yes
     63    2012-08-05 00:00:01  2012-08-05 00:00:02 ...    N/A      No
     64    2012-08-05 00:00:01  2012-08-05 00:00:02 ...    N/A      No
     65    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    N/A     Yes
     66    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    N/A     Yes
     67    2013-08-05 00:00:01  2013-08-05 00:00:02 ...    N/A      No
     68    2013-08-05 00:00:01  2013-08-05 00:00:02 ...    N/A      No
     69    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    N/A     Yes
     70    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    N/A     Yes
    Length = 38 rows

7. Querying the database
------------------------
The API for querying databases is similar to querying the VSO using the
method :meth:`sunpy.net.vso.VSOClient.search`. The :meth:`Database.search`
method accepts any number of ORed query attributes (using \|) and
combines them using AND. It returns a list of matched database entries.
The special thing about querying databases is that all attributes support
the unary operator ``~`` to negate specific attributes. Example: the query
``~Instrument('EIT')`` returns all entries that have *not* been observed
with the EIT.

7.1 Using VSO attributes
~~~~~~~~~~~~~~~~~~~~~~~~
Using the attributes from :mod:`sunpy.net.vso.attrs` is quite intuitive:
the simple attributes and the Time attribute work exactly as you expect.

.. note::

  The `near` parameter of :class:`sunpy.net.vso.attrs.Time` is ignored, because
  it's behaviour is not documented and it is different depending on the
  server which is requested.

The following query returns the data that was added in section 2.3.2:

    >>> print(display_entries(
    ...     database.search(vso.attrs.Time('2012-08-05', '2012-08-05 00:00:05'), vso.attrs.Instrument('AIA')),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax'], sort=True))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start observation_time_end instrument wavemin wavemax
    --- ---------------------- -------------------- ---------- ------- -------
     63    2012-08-05 00:00:01  2012-08-05 00:00:02        AIA     9.4     9.4
     64    2012-08-05 00:00:01  2012-08-05 00:00:02        AIA     9.4     9.4
     65    2012-08-05 00:00:02  2012-08-05 00:00:03        AIA    33.5    33.5
     66    2012-08-05 00:00:02  2012-08-05 00:00:03        AIA    33.5    33.5

.. NOTE the following code does not actually work. There seems to be a bug
    in sunpy.util.unit_conversion.to_angstrom (this is called by the
    __init__ of the Wave attribute)

When using the :class:`sunpy.net.vso.attrs.Wave` attribute, you have to
specify a unit using `astropy.units.Quantity`. If not an error is
raised. This also means that there is no default unit that is
used by the class. To know how you can specify a detail using astropy
check `astropy.units`.

    >>> from astropy import units as u
    >>> print(display_entries(
    ...     database.search(vso.attrs.Wavelength(1.0*u.nm, 2.0*u.nm)),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax'], sort=True))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start ...      wavemin            wavemax
    --- ---------------------- ... ------------------ ------------------
      1    2011-06-07 06:33:02 ...               17.1               17.1
     10    2011-06-07 06:39:31 ...               19.3               19.3
     11    2011-06-07 06:45:55 ...               19.3               19.3
     12    2011-06-07 06:52:19 ...               19.3               19.3
     13    2011-06-07 06:58:43 ...               19.3               19.3
     14    2011-06-07 20:37:52 ...               19.5               19.5
      2    2011-06-07 06:33:01 ... 13.100000000000001 13.100000000000001
     21    2011-06-07 06:33:29 ... 17.400000000000002 17.400000000000002
     22    2014-04-09 06:00:12 ...               17.1               17.1
      3    2011-06-07 06:33:02 ...               17.1               17.1
     59    2011-05-08 00:00:00 ...               17.1               17.1
      8    2011-06-07 06:33:07 ...               19.3               19.3
      9    2011-06-07 06:33:07 ...               19.3               19.3

7.2 Database-specific attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There are 5 additional query attributes supported by the database package.
They can be imported from the submodule :mod:`sunpy.database.attrs` and are in
particular:

- Starred
- Tag
- Path
- DownloadTime
- FitsHeaderEntry

The following query searches for all entries that have the tag 'spring' or
(inclusive or!) are starred and have not the FITS header key 'WAVEUNIT'
with the value 'Angstrom':

    >>> import sunpy.database.attrs as dbattrs
    >>> print(display_entries(
    ...     database.search(dbattrs.Tag('spring') | dbattrs.Starred(), ~dbattrs.FitsHeaderEntry('WAVEUNIT', 'Angstrom')),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax', 'tags', 'starred'], sort=True))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start observation_time_end ... wavemax  tags  starred
    --- ---------------------- -------------------- ... ------- ------ -------
     22    2014-04-09 06:00:12                  N/A ...    17.1 spring      No
      4    2011-06-07 06:33:02  2011-06-07 06:33:02 ...    21.1    N/A     Yes
      5    2011-06-07 06:33:03  2011-06-07 06:33:03 ...    33.5    N/A     Yes
     59    2011-05-08 00:00:00  2011-05-08 00:00:01 ...    17.1 spring      No
     60    2011-05-08 00:00:00  2011-05-08 00:00:01 ...    21.1 spring     Yes
     61    2011-05-08 00:00:02  2011-05-08 00:00:03 ...     9.4 spring      No
     62    2011-05-08 00:00:03  2011-05-08 00:00:04 ...    33.5 spring     Yes
     65    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    33.5    N/A     Yes
     66    2012-08-05 00:00:02  2012-08-05 00:00:03 ...    33.5    N/A     Yes
     69    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    33.5    N/A     Yes
      7    2011-06-07 06:33:05  2011-06-07 06:33:05 ...   160.0    N/A     Yes
     70    2013-08-05 00:00:02  2013-08-05 00:00:03 ...    33.5    N/A     Yes


8. Caching
----------
All entries that are saved in the database are also saved in a cache
in-memory. The type of the cache is determined at the initialization of
the database object and cannot be changed after that. The default type is
:class:`caching.LRUCache` (least-recently used) and the other one which is
supported is :class:`caching.LFUCache` (least-frequently used). By
default the cache size is ``float('inf')``, i.e. infinite. To set the
cache size after the database object has been initialized, use the method
:meth:`Database.set_cache_size`. If the new size is smaller than the
current number of database entries, entries are removed according to the
cache type until the number of entries is equal to the given cache size.
The following call to :meth:`Database.set_cache_size` sets the cache size
to 10 and therefore removes the 5 entries that been used least recently.

    >>> database.set_cache_size(10)
    >>> print(display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax']))   # doctest: +NORMALIZE_WHITESPACE +REMOTE_DATA
     id observation_time_start ...      wavemin            wavemax
     --- ---------------------- ... ------------------ ------------------
      8    2011-06-07 06:33:07 ...               19.3               19.3
      9    2011-06-07 06:33:07 ...               19.3               19.3
     10    2011-06-07 06:39:31 ...               19.3               19.3
     11    2011-06-07 06:45:55 ...               19.3               19.3
     12    2011-06-07 06:52:19 ...               19.3               19.3
     13    2011-06-07 06:58:43 ...               19.3               19.3
     14    2011-06-07 20:37:52 ...               19.5               19.5
     21    2011-06-07 06:33:29 ... 17.400000000000002 17.400000000000002
     22    2014-04-09 06:00:12 ...               17.1               17.1
     58    2011-06-06 00:00:00 ...                N/A                N/A
