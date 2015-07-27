.. _database_guide:

--------------------------
Using the database package
--------------------------
.. currentmodule:: sunpy.database

The database package offers the possibility to save retrieved data (e.g.
via the :mod:`sunpy.net.vso` package) onto a local or remote database. The
database may be a single file located on a local hard drive (if a SQLite
database is used) or a local or remote database server (see the SQLAlchemy
documentation for a `list of supported databases
<http://docs.sqlalchemy.org/en/rel_0_8/dialects/>`_)
This makes it possible to fetch required data from the local database
instead of downloading it again from a remote server. The package
:mod:`sunpy.database` was developed as part of `Google Summer of Code
(GSOC) 2013
<http://www.google-melange.com/gsoc/homepage/google/gsoc2013>`_.

1. Connecting and initializing the database
-------------------------------------------
To start a connection to an existing or a new database, instantiate
a :class:`Database` object.

    >>> from sunpy.database import Database
    >>> database = Database('sqlite:///sunpydata.sqlite')

The database object in our example above connects to a new SQLite database with
the file name "sunpydata.sqlite" in the current directory.

The first parameter of ``Database`` receives one mandatory argument:
a URL which describes how to connect to the database. This value is
directly passed to :func:`sqlalchemy.create_engine`. The supported
format of this URL is described by the documentation of
:func:`sqlalchemy.create_engine` as follows:

    "The string form of the URL is
    ``dialect+driver://user:password@host/dbname[?key=value..]``, where
    dialect is a database name such as ``mysql``, ``oracle``, ``postgresql``,
    etc., and driver the name of a DBAPI, such as ``psycopg2``,
    ``pyodbc``, ``cx_oracle``, etc."

Note that a connection is only established when it's really needed, i.e. if some query
is sent to the database to read from it. Transactions can also be committed explicitly
using the :meth:`Database.commit` method.

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
Before explaining how to add new entries, it is important to know what
information is saved in an entry. Each database entry is an instance of
the class :class:`tables.DatabaseEntry` with the following attributes:

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
  :meth:`Database.download` and also possibly by :meth:`Database.fetch`. ``starred``
  is set or changed via the method   :meth:`Database.star` or :meth:`unstar`,
  respectively. Analogously, ``tags`` is set via the methods :meth:`Database.tag`
  and :meth:`Database.remove_tag`.

* The attribute ``fits_header_entries`` is set by the methods
  :meth:`Database.download`, :meth:`Database.add_from_dir`, and
  :meth:`Database.add_from_file`.

2.1 Adding entries from one FITS file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The method :meth:`Database.add_from_file` receives one positional argument
(either a string or a file-like object) which is used to add at least one
new entry from the given FITS file to the database. Why "at least one" and
not "exactly one"? The reason is that each database entry does not
represent one file but one FITS header of a file. That means, if you pass
a file which has 5 FITS headers in it, 5 entries will be added to the
database. The file in the following example (``sunpy.data.sample.AIA_171_IMAGE``) has
only one FITS header, that is why just one entry is added to the database.
However, if you are working with Hinode/SOT files you may notice that for
each file you get two entries, one which refers to the observation and
another that contains some (useless - as `discussed
<http://listmgr.cv.nrao.edu/pipermail/fitsbits/2007-August/001923.html>`_
in the `fitsbits mailing list
<http://listmgr.cv.nrao.edu/mailman/listinfo/fitsbits>`_)
telemetry data.

The method saves the value of `path` by either simply passing on the value
of the received argument (if it was a string) or by reading the value of
``file.name`` where ``file`` is the passed argument. If the path cannot be
determined, it stays to ``None`` (the default value).

The values of `wavemin` and `wavemax` are only set if the wavelength unit
of the passed FITS file can be found out or if the attribute
`default_waveunit` of the database object is set. These values are then
used to convert from the used unit to nanometers. The rationale behind
this behaviour is that it makes querying for wavelengths more flexible.
tl;dr: **wavelengths are always stored in nanometers!**

The value of the attribute `instrument` is simply set by looking up the
FITS header key *INSTRUME*. The value of `observation_time_start` is set
by searching for the FITS header key *DATE-OBS* or *DATE_OBS*.
Analogously, `observation_time_end` is set by searching for *DATE-END* or
*DATE_END*. Finally, the whole FITS header is stored in the attribute
`fits_header_entries` as a list of :class:`tables.FitsHeaderEntry`
instances. All FITS comments are stored in the attribute
`fits_key_comments` which is a list of :class:`tables.FitsKeyComment`
instances.

Using the function ``len`` on a :class:`Database` object returns the
number of saved database entries. To get the first entry of the database,
``database[0]`` is used (the ID number of the entries does not matter,
``database[0]`` always returns the oldest saved entry of the database). If
the database had been empty, this expression would have raised an
:exc:`IndexError`. In section 3, more advanced formats of the slicing
syntax are introduced.

    >>> import sunpy.data
    >>> sunpy.data.download_sample_data(overwrite=False)  # doctest: +SKIP
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
    (datetime.datetime(2011, 3, 19, 10, 54, 0, 340000), None)
    >>> len(entry.fits_header_entries)
    170
    >>> for fits_header_entry in entry.fits_header_entries[:10]:
    ...     print '{entry.key}\n\t{entry.value}'.format(entry=fits_header_entry)   # doctest: +NORMALIZE_WHITESPACE
    SIMPLE
    	True
    BITPIX
    	32
    NAXIS
     	2
    NAXIS1
    	1024
    NAXIS2
        1024
    EXTEND
    	True
    COMMENT
    	FITS (Flexible Image Transport System) format is defined in 'Astronomy  and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
    ORIGIN
    	SDO/JSOC-SDP
    DATE
        2011-03-19T11:08:25
    TELESCOP
       	SDO/AIA

    >>> for fits_key_comment in entry.fits_key_comments:
    ...     print '{comment.key}\n\t{comment.value}'.format(comment=fits_key_comment)   # doctest: +NORMALIZE_WHITESPACE
    NAXIS
            number of data axes
    NAXIS1
            length of data axis 1
    DATASUM
            data unit checksum updated 2011-03-19T11:08:18
    EXTEND
            FITS dataset may contain extensions
    BITPIX
            number of bits per data pixel
    SIMPLE
            file does conform to FITS standard
    CHECKSUM
            HDU checksum updated 2011-03-19T11:08:18
    NAXIS2
            length of data axis 2


2.2 Adding entries from a directory of FITS files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Adding all FITS files from a certain directory works by calling the method
:meth:`Database.add_from_dir` and passing the desired directory to it. By
setting the keyword argument `ignore_already_added` to `True`, no
exception is raised if it is attempted to add an already existing entry
(In this case, setting this parameter is required because the file
``sunpy.data.sample.AIA_171_IMAGE`` was already added which is located in the
directory ``sampledata_dir``).

    >>> from sunpy import config
    >>> sampledata_dir = config.get("downloads", "sample_dir")
    >>> database.default_waveunit = 'angstrom'
    >>> database.add_from_dir(sampledata_dir, ignore_already_added=True)
    >>> len(database)
    25

2.3 Adding entries using the VSO interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

2.3.1 Adding entries from a VSO query result
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A VSO query result can be used to add new entries to the database. The
number of database entries that will be added is equal to the value of
``qr.num_records()`` in the following code snippet. Note that the method
:meth:`Database.add_from_vso_query_result` does not download any files,
though. If you want to add new entries using the VSO and also want to
download files at the same time, take a look at the following two
sections.

    >>> from sunpy.net import vso
    >>> client = vso.VSOClient()
    >>> qr = client.query(
    ...     vso.attrs.Time('2011-05-08', '2011-05-08 00:00:05'),
    ...     vso.attrs.Instrument('AIA'))
    >>> qr.num_records()
    4
    >>> database.add_from_vso_query_result(qr)
    >>> len(database)
    29

2.3.2 Downloading
^^^^^^^^^^^^^^^^^
The method :meth:`Database.download` queries the VSO by passing the given
query on to :meth:`sunpy.net.vso.VSOClient.query`. Note that not the
number of records of the resulting query result determines the number of
entries that will be added to the database! The number of entries that
will be added depends on the total number of FITS headers. The download
method also accepts an optional keyword argument `path` which is passed
as-is to :meth:`sunpy.net.vso.VSOClient.get` and determines the value of
the `path` attribute of each entry.

    >>> database.download(
    ...     vso.attrs.Time('2012-08-05', '2012-08-05 00:00:05'),
    ...     vso.attrs.Instrument('AIA'))
    >>> len(database)
    33

2.3.3 "Clever" Fetching
^^^^^^^^^^^^^^^^^^^^^^^
The method :meth:`Database.fetch` queries the database if the given query
has already been used once to add entries using the method
:meth:`Database.download`. Otherwise, the given query is used to download
and add new data via :meth:`Database.download`. This means: If you are not
sure whether the required data already exists in the database, use
:meth:`Database.fetch` and you use a method to save bandwidth!

As you can see in the following example, the first call of the fetch
method results in a list of entries. These are the entries which have been
returned by querying the database. This is the reason why the number of
saved entries in the database is still 32 even after the fetch method has
been called. The second fetch call has not been done already on a download
call, therefore the given query is used to download new data and add these
resulting entries to the database. Because the query result translates to
4 records, 4 new entries are added to the database after the fetch call.

    >>> entries = database.fetch(
    ...     vso.attrs.Time('2012-08-05', '2012-08-05 00:00:05'),
    ...     vso.attrs.Instrument('AIA'))
    >>> entries is None
    False
    >>> len(entries)
    4
    >>> len(database)
    33

    >>> entries = database.fetch(
    ...     vso.attrs.Time('2013-08-05', '2013-08-05 00:00:05'),
    ...     vso.attrs.Instrument('AIA'))
    >>> entries is None
    True
    >>> len(database)
    37

2.4 Adding entries manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Although usually not required, it is also possible to add database entries
by specifying the parameters manually. To do so, you simply pass the
values as keyword arguments to :class:`tables.DatabaseEntry` as follows:

    >>> from sunpy.database.tables import DatabaseEntry
    >>> entry = DatabaseEntry(instrument='EIT', wavemin=25.0)
    >>> database.add(entry)
    >>> entry in database
    False
    >>> database.commit()
    >>> entry in database
    True
    >>> len(database)
    38

Note that the `in` operator works only as expected after the
:meth:`Database.commit` method has been called!

3. Displaying entries in a table
--------------------------------
Meanwhile, 37 entries have been added, all of them saving a lot of data.
How can selected data of each entry be displayed? Fortunately, there is a
helper function to do this: :func:`tables.display_entries` takes two
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
    >>> print display_entries(database,
    ...                       ['id', 'observation_time_start', 'observation_time_end',
    ...                        'instrument', 'wavemin', 'wavemax'])   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax
    -- ---------------------- -------------------- ---------- ------- -------
    1  2011-03-19 10:54:00    N/A                  AIA_3      17.1    17.1
    2  N/A                    N/A                  N/A        N/A     N/A
    3  2013-09-21 16:00:06    N/A                  AIA_2      19.3    19.3
    4  2011-03-19 10:54:00    N/A                  AIA_3      17.1    17.1
    5  2014-04-09 06:00:12    N/A                  AIA_3      17.1    17.1
    6  2011-09-22 00:00:00    2011-09-22 00:00:00  BIR        N/A     N/A
    7  N/A                    N/A                  N/A        N/A     N/A
    8  2002-06-25 10:00:10    N/A                  EIT        19.5    19.5
    9  2002-02-20 11:06:00    2002-02-20 11:06:43  RHESSI     N/A     N/A
    10 N/A                    N/A                  N/A        N/A     N/A
    11 N/A                    N/A                  N/A        N/A     N/A
    12 N/A                    N/A                  N/A        N/A     N/A
    13 N/A                    N/A                  N/A        N/A     N/A
    14 N/A                    N/A                  N/A        N/A     N/A
    15 N/A                    N/A                  N/A        N/A     N/A
    16 N/A                    N/A                  N/A        N/A     N/A
    17 N/A                    N/A                  N/A        N/A     N/A
    18 N/A                    N/A                  N/A        N/A     N/A
    19 N/A                    N/A                  N/A        N/A     N/A
    20 2010-10-16 19:12:18    2010-10-16 19:12:22  RHESSI     N/A     N/A
    21 N/A                    N/A                  N/A        N/A     N/A
    22 N/A                    N/A                  N/A        N/A     N/A
    23 N/A                    N/A                  N/A        N/A     N/A
    24 2012-10-30 15:30:01    N/A                  AIA_4      9.4     9.4
    25 2012-01-01 00:16:07    N/A                  SWAP       17.4    17.4
    26 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        17.1    17.1
    27 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        21.1    21.1
    28 2011-05-08 00:00:02    2011-05-08 00:00:03  AIA        9.4     9.4
    29 2011-05-08 00:00:03    2011-05-08 00:00:04  AIA        33.5    33.5
    30 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4
    31 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4
    32 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5
    33 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5
    34 2013-08-05 00:00:01    2013-08-05 00:00:02  AIA        9.4     9.4
    35 2013-08-05 00:00:01    2013-08-05 00:00:02  AIA        9.4     9.4
    36 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5
    37 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5
    38 N/A                    N/A                  EIT        25.0    N/A

In Section 2.1, "Adding entries from one FITS file", it has already been
shows that the index operator can be used to access certain single
database entries like ``database[5]`` or even use negative indices such as
in ``database[-2]`` to get the second-latest entry. It is also possible to
use the more advanced slicing syntax: ``database[:]`` returns a list of
all database entries (and is hereby an alias for ``list(database)``),
``database[::2]`` returns a list of every 2nd entry, starting with the
oldest one. As you can imagine, ``database[9::10]`` starts with the 10th
entry and returns a list of every 10th entry from there.

    >>> print display_entries(database[9::10],
    ...                       ['id', 'observation_time_start', 'observation_time_end',
    ...                        'instrument', 'wavemin', 'wavemax'])   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax
    -- ---------------------- -------------------- ---------- ------- -------
    10 N/A                    N/A                  N/A        N/A     N/A
    20 2010-10-16 19:12:18    2010-10-16 19:12:22  RHESSI     N/A     N/A
    30 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4

4. Removing entries
-------------------
``database.remove()`` can be used to remove database entries from the SunPy
database. It takes a ``tables.DatabaseEntry`` object as argument.

For example, let us imagine we want to only have database entries which have some
observation time saved. To remove all entries where the value of both
``observation_time_start`` and ``observation_time_end`` is None, one can
simply iterate over the database and uses the :meth:`Database.remove`
method to remove those where the just described predicate is true:

    >>> for database_entry in database:
    ...     if database_entry.observation_time_start is None and database_entry.observation_time_end is None:
    ...         database.remove(database_entry)
    ...
    >>> len(database)
    22
    >>> print display_entries(database,
    ...                       ['id', 'observation_time_start', 'observation_time_end',
    ...                        'instrument', 'wavemin', 'wavemax'])   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax
    -- ---------------------- -------------------- ---------- ------- -------
    1  2011-03-19 10:54:00    N/A                  AIA_3      17.1    17.1
    3  2013-09-21 16:00:06    N/A                  AIA_2      19.3    19.3
    4  2011-03-19 10:54:00    N/A                  AIA_3      17.1    17.1
    5  2014-04-09 06:00:12    N/A                  AIA_3      17.1    17.1
    6  2011-09-22 00:00:00    2011-09-22 00:00:00  BIR        N/A     N/A
    8  2002-06-25 10:00:10    N/A                  EIT        19.5    19.5
    9  2002-02-20 11:06:00    2002-02-20 11:06:43  RHESSI     N/A     N/A
    20 2010-10-16 19:12:18    2010-10-16 19:12:22  RHESSI     N/A     N/A
    24 2012-10-30 15:30:01    N/A                  AIA_4      9.4     9.4
    25 2012-01-01 00:16:07    N/A                  SWAP       17.4    17.4
    26 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        17.1    17.1
    27 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        21.1    21.1
    28 2011-05-08 00:00:02    2011-05-08 00:00:03  AIA        9.4     9.4
    29 2011-05-08 00:00:03    2011-05-08 00:00:04  AIA        33.5    33.5
    30 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4
    31 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4
    32 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5
    33 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5
    34 2013-08-05 00:00:01    2013-08-05 00:00:02  AIA        9.4     9.4
    35 2013-08-05 00:00:01    2013-08-05 00:00:02  AIA        9.4     9.4
    36 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5
    37 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5

There are more possible ways to remove entries: You can remove every 2nd
entry by iterating over ``database[::2]`` and then calling passing every
yielded entry to :meth:`Database.remove`. You can even do advanced
operations easily like for example removing every entry with a certain
observation start time, instrument, and FITS header entry pair. This
requires knowledge of the :meth:`Database.query` method though, which will
be covered in section 7, "Querying the database".

5. Editing entries
------------------
5.1 Starring and unstarring entries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The database package supports marking certain entries as "starred". This
concept may be familiar to you from graphical applications such as E-Mail
clients or photo managers. The method :meth:`Database.star` marks the
passed entry as starred. Let's say we are for some reason interested in
all values that have a wavelength of 20nm or higher, so we mark those as
starred:

    >>> for database_entry in database:
    ...     if database_entry.wavemin > 20:
    ...         database.star(database_entry)
    >>> print display_entries(
    ...     filter(lambda entry: entry.starred, database),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax'])   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax
    -- ---------------------- -------------------- ---------- ------- -------
    27 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        21.1    21.1
    29 2011-05-08 00:00:03    2011-05-08 00:00:04  AIA        33.5    33.5
    32 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5
    33 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5
    36 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5
    37 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5


So remove the mark from these entries, the method :meth:`Database.unstar`
works the same way.

5.2 Setting and removing tags
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The starring concept is somewhat "binary" in a way that one can only
decide whether to mark an entry as starred or not. To add some more
information by assigning keywords to entries, the database package also
supports *tags*. The `tags` property of a database object holds all tags
that are saved in this database. Let us assign the tag *spring* to all
entries that have been observed in March, April, or May on any day at any
year:

    >>> for database_entry in database:
    ...     if database_entry.observation_time_start.month in [3,4,5]:
    ...         database.tag(database_entry, 'spring')
    >>> database.tags
    [<Tag(name 'spring')>]
    >>> spring = database.tags[0]
    >>> print display_entries(
    ...     filter(lambda entry: spring in entry.tags, database),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax'])   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax
    -- ---------------------- -------------------- ---------- ------- -------
    1  2011-03-19 10:54:00    N/A                  AIA_3      17.1    17.1
    4  2011-03-19 10:54:00    N/A                  AIA_3      17.1    17.1
    5  2014-04-09 06:00:12    N/A                  AIA_3      17.1    17.1
    26 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        17.1    17.1
    27 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        21.1    21.1
    28 2011-05-08 00:00:02    2011-05-08 00:00:03  AIA        9.4     9.4
    29 2011-05-08 00:00:03    2011-05-08 00:00:04  AIA        33.5    33.5

5.3 Changing custom attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
What sill annoys me a bit is that there are entries in the database where
the end of the observation time is `None`. Let us change them to the same
value to the start of the observation time (because we can and because it
looks prettier, not because it is accurate). The :meth:`Database.edit`
method receives the database entry to be edited and any number of keyword
arguments which describe which values to change and how.

    >>> for database_entry in database:
    ...     if database_entry.observation_time_end is None:
    ...         database.edit(database_entry, observation_time_end=database_entry.observation_time_start)
    ...
    >>> print display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax'])   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax
    -- ---------------------- -------------------- ---------- ------- -------
    1  2011-03-19 10:54:00    2011-03-19 10:54:00  AIA_3      17.1    17.1
    3  2013-09-21 16:00:06    2013-09-21 16:00:06  AIA_2      19.3    19.3
    4  2011-03-19 10:54:00    2011-03-19 10:54:00  AIA_3      17.1    17.1
    5  2014-04-09 06:00:12    2014-04-09 06:00:12  AIA_3      17.1    17.1
    6  2011-09-22 00:00:00    2011-09-22 00:00:00  BIR        N/A     N/A
    8  2002-06-25 10:00:10    2002-06-25 10:00:10  EIT        19.5    19.5
    9  2002-02-20 11:06:00    2002-02-20 11:06:43  RHESSI     N/A     N/A
    20 2010-10-16 19:12:18    2010-10-16 19:12:22  RHESSI     N/A     N/A
    24 2012-10-30 15:30:01    2012-10-30 15:30:01  AIA_4      9.4     9.4
    25 2012-01-01 00:16:07    2012-01-01 00:16:07  SWAP       17.4    17.4
    26 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        17.1    17.1
    27 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        21.1    21.1
    28 2011-05-08 00:00:02    2011-05-08 00:00:03  AIA        9.4     9.4
    29 2011-05-08 00:00:03    2011-05-08 00:00:04  AIA        33.5    33.5
    30 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4
    31 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4
    32 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5
    33 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5
    34 2013-08-05 00:00:01    2013-08-05 00:00:02  AIA        9.4     9.4
    35 2013-08-05 00:00:01    2013-08-05 00:00:02  AIA        9.4     9.4
    36 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5
    37 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5

You may ask yourself now "Why can't I simply use
``database_entry.observation_time_end = database_entry.observation_time_start``"?
Well, the answer is: you can, but it has one major disadvantage: you
cannot undo this operation if you don't use the methods of the database
objects such as :class:`Database.edit`. See the section 6, "Undoing and
redoing operations", to see how undoing and redoing works.

6. Undoing and redoing operations
---------------------------------
A very handy feature of the database package is that every operation can
be reverted that changes the database in some way. The Database class has
the methods :meth:`Database.undo()` and :meth:`Database.redo()` to undo
and redo the last n commands, respectively.

.. note::

    The undo and redo history are only saved in-memory! That means in
    particular, that if you work on a Database object in an interactive
    Python session and quit this session, the undo and redo history are
    lost.

In the following snippet, the operations from the sections 5.3, 5.2, and
5.1 are undone. Note the following changes: there are no longer any tags
anymore saved in the database, there is no entry which is starred and
there are again entries with no end of observation time.

    >>> database.undo(4)  # undo the edits from 5.3 (4 records have been affected)
    >>> database.undo(6)  # undo tagging some entries with the tag 'spring'
    >>> database.undo(4)  # undo starring entries
    >>> print display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax', 'tags', 'starred'])   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax tags starred
    -- ---------------------- -------------------- ---------- ------- ------- ---- -------
    1  2011-03-19 10:54:00    N/A                  AIA_3      17.1    17.1    N/A  No
    3  2013-09-21 16:00:06    N/A                  AIA_2      19.3    19.3    N/A  No
    4  2011-03-19 10:54:00    N/A                  AIA_3      17.1    17.1    N/A  No
    5  2014-04-09 06:00:12    N/A                  AIA_3      17.1    17.1    N/A  No
    6  2011-09-22 00:00:00    2011-09-22 00:00:00  BIR        N/A     N/A     N/A  No
    8  2002-06-25 10:00:10    N/A                  EIT        19.5    19.5    N/A  No
    9  2002-02-20 11:06:00    2002-02-20 11:06:43  RHESSI     N/A     N/A     N/A  No
    20 2010-10-16 19:12:18    2010-10-16 19:12:22  RHESSI     N/A     N/A     N/A  No
    24 2012-10-30 15:30:01    N/A                  AIA_4      9.4     9.4     N/A  No
    25 2012-01-01 00:16:07    N/A                  SWAP       17.4    17.4    N/A  No
    26 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        17.1    17.1    N/A  No
    27 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        21.1    21.1    N/A  Yes
    28 2011-05-08 00:00:02    2011-05-08 00:00:03  AIA        9.4     9.4     N/A  No
    29 2011-05-08 00:00:03    2011-05-08 00:00:04  AIA        33.5    33.5    N/A  Yes
    30 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4     N/A  No
    31 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4     N/A  No
    32 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5    N/A  Yes
    33 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5    N/A  Yes
    34 2013-08-05 00:00:01    2013-08-05 00:00:02  AIA        9.4     9.4     N/A  No
    35 2013-08-05 00:00:01    2013-08-05 00:00:02  AIA        9.4     9.4     N/A  No
    36 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5    N/A  Yes
    37 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5    N/A  Yes


The redo method reverts the last n operations that have been undone. If
not that many operations can be redone (i.e. any number greater than 14 in
this example), an exception is raised. You can see that the redo call
reverts the original state: the tags appeared again and all entries with a
wavelength >20nm are starred again. Also, there are no entries with no
stored end of observation time anymore.

    >>> database.redo(14)  # redo all undone operations
    >>> print display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax', 'tags', 'starred'])   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax tags   starred
    -- ---------------------- -------------------- ---------- ------- ------- ----   -------
    1  2011-03-19 10:54:00    2011-03-19 10:54:00  AIA_3      17.1    17.1    spring No
    3  2013-09-21 16:00:06    2013-09-21 16:00:06  AIA_2      19.3    19.3    N/A    No
    4  2011-03-19 10:54:00    2011-03-19 10:54:00  AIA_3      17.1    17.1    spring No
    5  2014-04-09 06:00:12    2014-04-09 06:00:12  AIA_3      17.1    17.1    spring No
    6  2011-09-22 00:00:00    2011-09-22 00:00:00  BIR        N/A     N/A     N/A    No
    8  2002-06-25 10:00:10    2002-06-25 10:00:10  EIT        19.5    19.5    N/A    No
    9  2002-02-20 11:06:00    2002-02-20 11:06:43  RHESSI     N/A     N/A     N/A    No
    20 2010-10-16 19:12:18    2010-10-16 19:12:22  RHESSI     N/A     N/A     N/A    No
    24 2012-10-30 15:30:01    2012-10-30 15:30:01  AIA_4      9.4     9.4     N/A    No
    25 2012-01-01 00:16:07    2012-01-01 00:16:07  SWAP       17.4    17.4    N/A    No
    26 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        17.1    17.1    spring No
    27 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        21.1    21.1    spring Yes
    28 2011-05-08 00:00:02    2011-05-08 00:00:03  AIA        9.4     9.4     spring No
    29 2011-05-08 00:00:03    2011-05-08 00:00:04  AIA        33.5    33.5    spring Yes
    30 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4     N/A    No
    31 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4     N/A    No
    32 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5    N/A    Yes
    33 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5    N/A    Yes
    34 2013-08-05 00:00:01    2013-08-05 00:00:02  AIA        9.4     9.4     N/A    No
    35 2013-08-05 00:00:01    2013-08-05 00:00:02  AIA        9.4     9.4     N/A    No
    36 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5    N/A    Yes
    37 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5    N/A    Yes

7. Querying the database
------------------------
The API for querying databases is similar to querying the VSO using the
method :meth:`sunpy.net.vso.VSOClient.query`. The :meth:`Database.query`
method accepts any number of ORed query attributes (using \|) and
combines them using AND. It returns a list of matched database entries.
The special thing about querying databases is that all attributes support
the unary operator ``~`` to negate specific attributes. Example: the query
``~Instrument('EIT')`` returns all entries that have *not* been observed
with the EIT.

7.1 Using VSO attributes
~~~~~~~~~~~~~~~~~~~~~~~~
Using the attributes from :mod:`sunpy.net.vso.attrs` is quite intuitive:
the simple attributes and the Time attribute work exactly as you expect
it. Note though that the `near` parameter of
:class:`sunpy.net.vso.attrs.Time` is ignored! The reason for this is that
its behaviour is not documented and that it is different depending on the
server which is requested. The following query returns the data that was
added in section 2.3.2, "Downloading":

    >>> print display_entries(
    ...     database.query(vso.attrs.Time('2012-08-05', '2012-08-05 00:00:05'), vso.attrs.Instrument('AIA')),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax'], sort=True)   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax
    -- ---------------------- -------------------- ---------- ------- -------
    30 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4
    31 2012-08-05 00:00:01    2012-08-05 00:00:02  AIA        9.4     9.4
    32 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5
    33 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5

.. NOTE the following code does not actually work. There seems to be a bug
    in sunpy.util.unit_conversion.to_angstrom (this is called by the
    __init__ of the Wave attribute)

When using the :class:`sunpy.net.vso.attrs.Wave` attribute, you have to
specify a unit using ``astropy.units.Quantity``. If not an error is
raised. This also implies that there is no default unit that is
used by the class. To know how you can specify a detail using astropy
check `astropy.units <https://astropy.readthedocs.org/en/stable/units/index.html>`_.

    >>> from astropy import units as u
    >>> print display_entries(
    ...     database.query(vso.attrs.Wave(1.0*u.nm, 2.0*u.nm)),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax'], sort=True)   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax
    -- ---------------------- -------------------- ---------- ------- -------
    1  2011-03-19 10:54:00    2011-03-19 10:54:00  AIA_3      17.1    17.1
    25 2012-01-01 00:16:07    2012-01-01 00:16:07  SWAP       17.4    17.4
    26 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        17.1    17.1
    3  2013-09-21 16:00:06    2013-09-21 16:00:06  AIA_2      19.3    19.3
    4  2011-03-19 10:54:00    2011-03-19 10:54:00  AIA_3      17.1    17.1
    5  2014-04-09 06:00:12    2014-04-09 06:00:12  AIA_3      17.1    17.1
    8  2002-06-25 10:00:10    2002-06-25 10:00:10  EIT        19.5    19.5


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
    >>> print display_entries(
    ...     database.query(dbattrs.Tag('spring') | dbattrs.Starred(), ~dbattrs.FitsHeaderEntry('WAVEUNIT', 'Angstrom')),
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax', 'tags', 'starred'], sort=True)   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax tags   starred
    -- ---------------------- -------------------- ---------- ------- ------- ----   -------
    1  2011-03-19 10:54:00    2011-03-19 10:54:00  AIA_3      17.1    17.1    spring No
    26 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        17.1    17.1    spring No
    27 2011-05-08 00:00:00    2011-05-08 00:00:01  AIA        21.1    21.1    spring Yes
    28 2011-05-08 00:00:02    2011-05-08 00:00:03  AIA        9.4     9.4     spring No
    29 2011-05-08 00:00:03    2011-05-08 00:00:04  AIA        33.5    33.5    spring Yes
    32 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5    N/A    Yes
    33 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5    N/A    Yes
    36 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5    N/A    Yes
    37 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5    N/A    Yes
    4  2011-03-19 10:54:00    2011-03-19 10:54:00  AIA_3      17.1    17.1    spring No
    5  2014-04-09 06:00:12    2014-04-09 06:00:12  AIA_3      17.1    17.1    spring No


8. Caching
----------
All entries that are saved in the database are also saved in a cache
in-memory. The type of the cache is determined at the initialization of
the database object and cannot be changed after that. The default type is
:class:`caching.LRUCache` (least-recently used) and the other one which is
supported is :class:`caching.LFUCache` (least-frequently used). Per
default the cache size is ``float('inf')``, i.e. infinite. To set the
cache size after the database object has been initialized, use the method
:meth:`Database.set_cache_size`. If the new size is smaller than the
current number of database entries, entries are removed according to the
cache type until the number of entries is equal to the given cache size.
The following call to :meth:`Database.set_cache_size` sets the cache size
to 10 and therefore removes the 5 entries that been used least recently.

    >>> database.set_cache_size(10)
    >>> print display_entries(
    ...     database,
    ...     ['id', 'observation_time_start', 'observation_time_end',
    ...      'instrument', 'wavemin', 'wavemax'])   # doctest: +NORMALIZE_WHITESPACE +SKIP
    id observation_time_start observation_time_end instrument wavemin wavemax
    -- ---------------------- -------------------- ---------- ------- -------
    1  2011-03-19 10:54:00    2011-03-19 10:54:00  AIA_3      17.1    17.1
    3  2013-09-21 16:00:06    2013-09-21 16:00:06  AIA_2      19.3    19.3
    4  2011-03-19 10:54:00    2011-03-19 10:54:00  AIA_3      17.1    17.1
    5  2014-04-09 06:00:12    2014-04-09 06:00:12  AIA_3      17.1    17.1
    8  2002-06-25 10:00:10    2002-06-25 10:00:10  EIT        19.5    19.5
    24 2012-10-30 15:30:01    2012-10-30 15:30:01  AIA_4      9.4     9.4
    25 2012-01-01 00:16:07    2012-01-01 00:16:07  SWAP       17.4    17.4
    33 2012-08-05 00:00:02    2012-08-05 00:00:03  AIA        33.5    33.5
    36 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5
    37 2013-08-05 00:00:02    2013-08-05 00:00:03  AIA        33.5    33.5
