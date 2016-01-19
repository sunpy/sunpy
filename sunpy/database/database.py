# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

from __future__ import absolute_import

import itertools
import operator
from datetime import datetime
from contextlib import contextmanager
import os.path

from sqlalchemy import create_engine, exists
from sqlalchemy.orm import sessionmaker

import sunpy
from sunpy.database import commands, tables, serialize
from sunpy.database.caching import LRUCache
from sunpy.database.commands import CompositeOperation
from sunpy.database.attrs import walker
from sunpy.net.hek2vso import H2VClient
from sunpy.net.attr import and_
from sunpy.net.vso import VSOClient

__authors__ = ['Simon Liedtke', 'Rajul Srivastava']
__emails__ = [
    'liedtke.simon@googlemail.com',
    'rajul09@gmail.com'
]


class EntryNotFoundError(Exception):
    """This exception is raised if a database entry cannot be found by its
    unique ID.

    """
    def __init__(self, entry_id):
        self.entry_id = entry_id

    def __str__(self):  # pragma: no cover
        return 'an entry with the ID {0:d} does not exist'.format(
            self.entry_id)


class EntryAlreadyAddedError(Exception):
    """This exception is raised if a database entry is attempted to be added to
    the database although it was already saved in it.

    """
    def __init__(self, database_entry):
        self.database_entry = database_entry

    def __str__(self):  # pragma: no cover
        return (
            'the entry {0!r} was already added '
            'to the database'.format(self.database_entry))


class EntryAlreadyStarredError(Exception):
    """This exception is raised if a database entry is marked as starred
    using :meth:`Database.star` although it was already starred before this
    operation.

    """
    def __init__(self, database_entry):
        self.database_entry = database_entry

    def __str__(self):  # pragma: no cover
        return (
            'the entry {0!r} is already marked '
            'as starred'.format(self.database_entry))


class EntryAlreadyUnstarredError(Exception):
    """This exception is raised if the star mark from a database entry is
    attempted to be removed although the entry is not starred.

    """
    def __init__(self, database_entry):
        self.database_entry = database_entry

    def __str__(self):  # pragma: no cover
        return (
            'the entry {0!r} is already not marked '
            'as starred'.format(self.database_entry))


class NoSuchTagError(Exception):
    """This exception is raised if a tag cannot be found in a database by its
    name.

    """
    def __init__(self, tag_name):
        self.tag_name = tag_name

    def __str__(self):  # pragma: no cover
        return 'the tag {0!r} is not saved in the database'.format(
            self.tag_name)


class TagAlreadyAssignedError(Exception):
    """This exception is raised if it is attempted to assign a tag to a
    database entry but the database entry already has this tag assigned.

    """
    def __init__(self, database_entry, tag_name):
        self.database_entry = database_entry
        self.tag_name = tag_name

    def __str__(self):  # pragma: no cover
        errmsg = 'the database entry {0!r} has already assigned the tag {1!r}'
        return errmsg.format(self.database_entry, self.tag_name)


@contextmanager
def disable_undo(database):
    """A context manager to disable saving the used commands in the undo
    history. This may be useful when it's important to save memory because a
    big number of entries in the undo history may occupy a lot of memory space.

    Examples
    --------
    >>> from sunpy.database import disable_undo, Database
    >>> from sunpy.database.tables import DatabaseEntry
    >>> database = Database('sqlite:///:memory:')
    >>> entry = DatabaseEntry()
    >>> with disable_undo(database) as db:
    ...     db.add(entry)

    # This will raise an EmptyCommandStackError
    >>> database.undo()   # doctest: +SKIP
    """
    database._enable_history = False
    yield database
    database._enable_history = True


class Database(object):
    """
    Database(url[, CacheClass[, cache_size[, default_waveunit]]])

    Parameters
    ----------
    url : str
        A URL describing the database. This value is simply passed to
        :func:`sqlalchemy.create_engine`
        If not specified the value will be read from the sunpy config file.
    CacheClass : sunpy.database.caching.BaseCache
        A concrete cache implementation of the abstract class BaseCache.
        Builtin supported values for this parameters are
        :class:`sunpy.database.caching.LRUCache` and
        :class:`sunpy.database.caching.LFUCache`.
        The default value is :class:`sunpy.database.caching.LRUCache`.
    cache_size : int
        The maximum number of database entries, default is no limit.
    default_waveunit : str, optional
        The wavelength unit that will be used if an entry is added to the
        database but its wavelength unit cannot be found (either in the file or
        the VSO query result block, depending on the way the entry was added).
        If `None` (the default), attempting to add an entry without knowing the
        wavelength unit results in a
        :exc:`sunpy.database.WaveunitNotFoundError`.
    """
    """
    Attributes
    ----------
    session : sqlalchemy.orm.session.Session
        A SQLAlchemy session object. This may be used for advanced queries and
        advanced manipulations and should only be used by people who are
        experienced with SQLAlchemy.

    cache_size: int
        The maximum number of database entries. This attribute is read-only. To
        change this value, use the method
        :meth:`sunpy.database.Database.set_cache_size`.

    tags : list of sunpy.database.Tag objects
        A list of all saved tags in database. This attribute is read-only.

    default_waveunit : str
        See "Parameters" section.

    Methods
    -------
    set_cache_size(cache_size)
        Set a new value for the maximum number of database entries in the
        cache. Use the value ``float('inf')`` to disable caching.
    commit()
        Flush pending changes and commit the current transaction.
    get_entry_by_id(id)
        Get the database entry which has the given unique ID number assigned.
    get_tag(tagname)
        Get the tag which has the given unique tagname assigned. Returns None
        if no tag with the given name is saved in the database.
    tag(entry, *tags)
        Assign the given database entry the given tags. If no tags are given,
        TypeError is raised.
    star(entry, ignore_already_starred=False)
        Mark the given database entry as starred. If ``ignore_already_starred``
        is False and the given entry is already marked as starred,
        EntryAlreadyStarredError is raised.
    unstar(entry, ignore_already_unstarred=False)
        Remove the starred mark of the given entry. If
        ``ignore_already_unstarred`` is False and the entry is not marked as
        starred, EntryAlreadyUnstarredError is raised.
    add(entry, ignore_already_added=False)
        Add the given database entry to the database. If
        ``ignore_already_added`` is False and the given entry is already saved
        in the database, EntryAlreadyAddedError is raised.
    edit(entry, **kwargs)
        Change the given database entry so that it interprets the passed
        key-value pairs as new values where the keys represent the attributes
        of this entry. If no keywords arguments are given, :exc:`ValueError` is
        raised.
    remove(entry)
        Remove the given entry from the database.
    undo(n=1)
        Redo the last n operations.
    redo(n=1)
        Redo the last n undone operations.
    __contains__(entry)
        Return True if the given database entry is saved in the database,
        False otherwise.
    __iter__()
        Return an iterator over all database entries.
    __len__()
        Get the number of database entries.

    """
    def __init__(self, url=None, CacheClass=LRUCache, cache_size=float('inf'),
            default_waveunit=None):
        if url is None:
            url = sunpy.config.get('database', 'url')
        self._engine = create_engine(url)
        self._session_cls = sessionmaker(bind=self._engine)
        self.session = self._session_cls()
        self._command_manager = commands.CommandManager()
        self.default_waveunit = default_waveunit
        self._enable_history = True

        class Cache(CacheClass):
            def callback(this, entry_id, database_entry):
                self.remove(database_entry)

            def append(this, value):
                this[max(this or [0]) + 1] = value
        self._create_tables()
        self._cache = Cache(cache_size)
        for entry in self:
            self._cache[entry.id] = entry

    @property
    def url(self):
        """The sqlalchemy url of the database instance"""
        return str(self._engine.url)

    @property
    def cache_size(self):
        return len(self._cache)

    @property
    def cache_maxsize(self):
        return self._cache.maxsize

    def set_cache_size(self, cache_size):
        """Set a new value for the maximum number of database entries in the
        cache. Use the value ``float('inf')`` to disable caching. If the new
        cache is smaller than the previous one and cannot contain all the
        entries anymore, entries are removed from the cache until the number of
        entries equals the cache size. Which entries are removed depends on the
        implementation of the cache (e.g.
        :class:`sunpy.database.caching.LRUCache`,
        :class:`sunpy.database.caching.LFUCache`).

        """
        cmds = CompositeOperation()
        # remove items from the cache if the given argument is lower than the
        # current cache size
        while cache_size < self.cache_size:
            # remove items from the cache until cache_size == maxsize of the
            # cache
            entry_id, entry = self._cache.to_be_removed
            cmd = commands.RemoveEntry(self.session, entry)
            if self._enable_history:
                cmds.add(cmd)
            else:
                cmd()
            del self._cache[entry_id]
        self._cache.maxsize = cache_size
        if cmds:
            self._command_manager.do(cmds)

    def _create_tables(self, checkfirst=True):
        """Initialise the database by creating all necessary tables. If
        ``checkfirst`` is True, already existing tables are not attempted to be
        created.

        """
        metadata = tables.Base.metadata
        metadata.create_all(self._engine, checkfirst=checkfirst)

    def commit(self):
        """Flush pending changes and commit the current transaction. This is a
        shortcut for :meth:`session.commit()`.

        """
        self.session.commit()

    def _download_and_collect_entries(self, query_result, client=None,
            path=None, progress=False):
        if client is None:
            client = VSOClient()

        paths = client.get(query_result, path).wait(progress=progress)

        for (path, block) in zip(paths, query_result):
            qr_entry = tables.DatabaseEntry._from_query_result_block(block)

            if os.path.isfile(path):
                entries = tables.entries_from_file(path, self.default_waveunit)
            elif os.path.isdir(path):
                entries = tables.entries_from_dir(path, self.default_waveunit)
            else:
                raise ValueError('The path is neither a file nor directory')

            for entry in entries:
                entry.source = qr_entry.source
                entry.provider = qr_entry.provider
                entry.physobs = qr_entry.physobs
                entry.fileid = qr_entry.fileid
                entry.observation_time_start = qr_entry.observation_time_start
                entry.observation_time_end = qr_entry.observation_time_end
                entry.instrument = qr_entry.instrument
                entry.size = qr_entry.size
                entry.wavemin = qr_entry.wavemin
                entry.wavemax = qr_entry.wavemax
                entry.path = path
                entry.download_time = datetime.utcnow()
                yield entry

    def download(self, *query, **kwargs):
        """download(*query, client=sunpy.net.vso.VSOClient(), path=None, progress=False)
        Search for data using the VSO interface (see
        :meth:`sunpy.net.vso.VSOClient.query`). If querying the VSO results in
        no data, no operation is performed. Concrete, this means that no entry
        is added to the database and no file is downloaded. Otherwise, the
        retrieved search result is used to download all files that belong to
        this search result. After that, all the gathered information (the one
        from the VSO query result and the one from the downloaded FITS files)
        is added to the database in a way that each FITS header is represented
        by one database entry.

        """
        if not query:
            raise TypeError('at least one attribute required')
        client = kwargs.pop('client', None)
        path = kwargs.pop('path', None)
        progress = kwargs.pop('progress', False)
        if kwargs:
            k, v = kwargs.popitem()
            raise TypeError('unexpected keyword argument {0!r}'.format(k))
        if client is None:
            client = VSOClient()
        qr = client.query(*query)
        # don't do anything if querying the VSO results in no data
        if not qr:
            return
        entries = list(self._download_and_collect_entries(
            qr, client, path, progress))
        dump = serialize.dump_query(and_(*query))
        (dump_exists,), = self.session.query(
            exists().where(tables.JSONDump.dump == tables.JSONDump(dump).dump))
        if dump_exists:
            # dump already exists in table jsondumps -> edit instead of add
            # update all entries with the fileid `entry.fileid`
            for entry in entries:
                old_entry = self.session.query(
                    tables.DatabaseEntry).filter_by(fileid=entry.fileid).first()
                if old_entry is not None:
                    attrs = [
                        'source', 'provider', 'physobs',
                        'observation_time_start', 'observation_time_end',
                        'instrument', 'size', 'wavemin', 'wavemax',
                        'download_time']
                    kwargs = dict((k, getattr(entry, k)) for k in attrs)
                    cmd = commands.EditEntry(old_entry, **kwargs)
                    if self._enable_history:
                        self._command_manager.do(cmd)
                    else:
                        cmd()
        else:
            self.add_many(entries)
            # serialize the query and save the serialization in the database
            # for two reasons:
            #   1. to avoid unnecessary downloading in future calls of
            #      ``fetch``
            #   2. to know whether to add or to edit entries in future calls of
            #      ``download`` (this method)
            self.session.add(tables.JSONDump(dump))

    def fetch(self, *query, **kwargs):
        """fetch(*query[, path])
        Check if the query has already been used to collect new data using the
        :meth:`sunpy.database.Database.download` method. If yes, query the
        database using the method :meth:`sunpy.database.Database.query` and
        return the result. Otherwise, call
        :meth:`sunpy.database.Database.download` and return the result.

        """
        if not query:
            raise TypeError('at least one attribute required')
        path = kwargs.pop('path', None)
        if kwargs:
            k, v = kwargs.popitem()
            raise TypeError('unexpected keyword argument {0!r}'.format(k))
        dump = serialize.dump_query(and_(*query))
        (dump_exists,), = self.session.query(
            exists().where(tables.JSONDump.dump == tables.JSONDump(dump).dump))
        if dump_exists:
            return self.query(*query)
        return self.download(*query, path=path)

    def query(self, *query, **kwargs):
        """
        query(*query[, sortby])
        Send the given query to the database and return a list of
        database entries that satisfy all of the given attributes.

        Apart from the attributes supported by the VSO interface, the following
        attributes are supported:

            - :class:`sunpy.database.attrs.Starred`

            - :class:`sunpy.database.attrs.Tag`

            - :class:`sunpy.database.attrs.Path`

            - :class:`sunpy.database.attrs.DownloadTime`

            - :class:`sunpy.database.attrs.FitsHeaderEntry`

        An important difference to the VSO attributes is that these attributes
        may also be used in negated form using the tilde ~ operator.

        Parameters
        ----------
        query : list
            A variable number of attributes that are chained together via the
            boolean AND operator. The | operator may be used between attributes
            to express the boolean OR operator.
        sortby : str, optional
            The column by which to sort the returned entries. The default is to
            sort by the start of the observation. See the attributes of
            :class:`sunpy.database.tables.DatabaseEntry` for a list of all
            possible values.

        Raises
        ------
        TypeError
            if no attribute is given or if some keyword argument other than
            'sortby' is given.

        Examples
        --------
        The query in the following example searches for all non-starred entries
        with the tag 'foo' or 'bar' (or both).

        >>> database.query(~attrs.Starred(), attrs.Tag('foo') | attrs.Tag('bar'))   # doctest: +SKIP

        """
        if not query:
            raise TypeError('at least one attribute required')
        sortby = kwargs.pop('sortby', 'observation_time_start')
        if kwargs:
            k, v = kwargs.popitem()
            raise TypeError('unexpected keyword argument {0!r}'.format(k))
        return sorted(
            walker.create(and_(*query), self.session),
            key=operator.attrgetter(sortby))

    def get_entry_by_id(self, entry_id):
        """Get a database entry by its unique ID number. If an entry with the
        given ID does not exist, :exc:`sunpy.database.EntryNotFoundError` is
        raised.

        """
        try:
            return self._cache[entry_id]
        except KeyError:
            raise EntryNotFoundError(entry_id)

    @property
    def tags(self):
        return self.session.query(tables.Tag).all()

    def get_tag(self, tag_name):
        """Get the tag which has the given name. If no such tag exists,
        :exc:`sunpy.database.NoSuchTagError` is raised.

        """
        for tag in self.tags:
            if tag_name == tag.name:
                return tag
        raise NoSuchTagError(tag_name)

    def tag(self, database_entry, *tags):
        """Assign the given database entry the given tags.

        Raises
        ------
        TypeError
            If no tags are given.

        sunpy.database.TagAlreadyAssignedError
            If at least one of the given tags is already assigned to the given
            database entry.

        """
        if not tags:
            raise TypeError('at least one tag must be given')
        # avoid duplicates
        tag_names = set(tags)
        cmds = CompositeOperation()
        for tag_name in tag_names:
            try:
                tag = self.get_tag(tag_name)
                if tag in database_entry.tags:
                    raise TagAlreadyAssignedError(database_entry, tag_names)
            except NoSuchTagError:
                # tag does not exist yet -> create it
                tag = tables.Tag(tag_name)
            cmd = commands.AddTag(self.session, database_entry, tag)
            if self._enable_history:
                cmds.add(cmd)
            else:
                cmd()
        if cmds:
            self._command_manager.do(cmds)

    def remove_tag(self, database_entry, tag_name):
        """Remove the given tag from the database entry. If the tag is not
        connected to any entry after this operation, the tag itself is removed
        from the database as well.

        Raises
        ------
        sunpy.database.NoSuchTagError
            If the tag is not connected to the given entry.

        """
        tag = self.get_tag(tag_name)
        cmds = CompositeOperation()
        remove_tag_cmd = commands.RemoveTag(self.session, database_entry, tag)
        remove_tag_cmd()
        if self._enable_history:
            cmds.add(remove_tag_cmd)
        if not tag.data:
            remove_entry_cmd = commands.RemoveEntry(self.session, tag)
            remove_entry_cmd()
            if self._enable_history:
                cmds.add(remove_entry_cmd)
        if self._enable_history:
            self._command_manager.push_undo_command(cmds)

    def star(self, database_entry, ignore_already_starred=False):
        """Mark the given database entry as starred. If this entry is already
        marked as starred, the behaviour depends on the optional argument
        ``ignore_already_starred``: if it is ``False`` (the default),
        :exc:`sunpy.database.EntryAlreadyStarredError` is raised. Otherwise,
        the entry is kept as starred and no exception is raised.

        """
        if database_entry.starred and not ignore_already_starred:
            raise EntryAlreadyStarredError(database_entry)
        self.edit(database_entry, starred=True)

    def unstar(self, database_entry, ignore_already_unstarred=False):
        """Remove the starred mark of the given entry. If this entry is not
        marked as starred, the behaviour depends on the optional argument
        ``ignore_already_unstarred``: if it is ``False`` (the default),
        :exc:`sunpy.database.EntryAlreadyUnstarredError` is raised. Otherwise,
        the entry is kept as unstarred and no exception is raised.

        """
        if not database_entry.starred and not ignore_already_unstarred:
            raise EntryAlreadyUnstarredError(database_entry)
        self.edit(database_entry, starred=False)

    def add_many(self, database_entries, ignore_already_added=False):
        """Add a row of database entries "at once". If this method is used,
        only one entry is saved in the undo history.

        Parameters
        ----------
        database_entries : iterable of sunpy.database.tables.DatabaseEntry
            The database entries that will be added to the database.

        ignore_already_added : bool, optional
            See Database.add

        """
        cmds = CompositeOperation()
        for database_entry in database_entries:
            # use list(self) instead of simply self because __contains__ checks
            # for existence in the database and not only all attributes except
            # ID.
            if database_entry in list(self) and not ignore_already_added:
                raise EntryAlreadyAddedError(database_entry)
            cmd = commands.AddEntry(self.session, database_entry)
            if self._enable_history:
                cmds.add(cmd)
            else:
                cmd()
            if database_entry.id is None:
                self._cache.append(database_entry)
            else:
                self._cache[database_entry.id] = database_entry
        if cmds:
            self._command_manager.do(cmds)

    def add(self, database_entry, ignore_already_added=False):
        """Add the given database entry to the database table.

        Parameters
        ----------
        database_entry : sunpy.database.tables.DatabaseEntry
            The database entry that will be added to this database.

        ignore_already_added : bool, optional
            If True, attempts to add an already existing database entry will
            result in a :exc:`sunpy.database.EntryAlreadyAddedError`.
            Otherwise, a new entry will be added and there will be duplicates
            in the database.

        """
        if database_entry in self and not ignore_already_added:
            raise EntryAlreadyAddedError(database_entry)
        add_entry_cmd = commands.AddEntry(self.session, database_entry)
        if self._enable_history:
            self._command_manager.do(add_entry_cmd)
        else:
            add_entry_cmd()
        if database_entry.id is None:
            self._cache.append(database_entry)
        else:
            self._cache[database_entry.id] = database_entry

    def add_from_hek_query_result(self, query_result,
            ignore_already_added=False):
        """Add database entries from a HEK query result.

        Parameters
        ----------
        query_result : list
            The value returned by :meth:`sunpy.net.hek.HEKClient().query`

        ignore_already_added : bool
            See :meth:`sunpy.database.Database.add`.

        """
        vso_qr = itertools.chain.from_iterable(
            H2VClient().translate_and_query(query_result))
        self.add_from_vso_query_result(vso_qr, ignore_already_added)

    def download_from_vso_query_result(self, query_result, client=None,
            path=None, progress=False, ignore_already_added=False):
        """download(query_result, client=sunpy.net.vso.VSOClient(),
        path=None, progress=False, ignore_already_added=False)

        Add new database entries from a VSO query result and download the
        corresponding data files. See :meth:`sunpy.database.Database.download`
        for information about the parameters `client`, `path`, `progress`.

        Parameters
        ----------
        query_result : sunpy.net.vso.vso.QueryResponse
            A VSO query response that was returned by the ``query`` method of a
            :class:`sunpy.net.vso.VSOClient` object.

        ignore_already_added : bool
            See :meth:`sunpy.database.Database.add`.

        """
        if not query_result:
            return
        self.add_many(self._download_and_collect_entries(
            query_result, client, path, progress))

    def add_from_vso_query_result(self, query_result,
            ignore_already_added=False):
        """Generate database entries from a VSO query result and add all the
        generated entries to this database.

        Parameters
        ----------
        query_result : sunpy.net.vso.vso.QueryResponse
            A VSO query response that was returned by the ``query`` method of a
            :class:`sunpy.net.vso.VSOClient` object.

        ignore_already_added : bool
            See :meth:`sunpy.database.Database.add`.

        """
        self.add_many(
            tables.entries_from_query_result(
                query_result, self.default_waveunit),
            ignore_already_added)

    def add_from_dir(self, path, recursive=False, pattern='*',
            ignore_already_added=False):
        """Search the given directory for FITS files and use their FITS headers
        to add new entries to the database. Note that one entry in the database
        is assigned to a list of FITS headers, so not the number of FITS headers
        but the number of FITS files which have been read determine the number
        of database entries that will be added. FITS files are detected by
        reading the content of each file, the `pattern` argument may be used to
        avoid reading entire directories if one knows that all FITS files have
        the same filename extension.

        Parameters
        ----------
        path : string
            The directory where to look for FITS files.

        recursive : bool, optional
            If True, the given directory will be searched recursively.
            Otherwise, only the given directory and no subdirectories are
            searched. The default is `False`, i.e. the given directory is not
            searched recursively.

        pattern : string, optional
            The pattern can be used to filter the list of filenames before the
            files are attempted to be read. The default is to collect all
            files. This value is passed to the function :func:`fnmatch.filter`,
            see its documentation for more information on the supported syntax.

        ignore_already_added : bool, optional
            See :meth:`sunpy.database.Database.add`.

        """
        cmds = CompositeOperation()
        entries = tables.entries_from_dir(
            path, recursive, pattern, self.default_waveunit)
        for database_entry, filepath in entries:
            if database_entry in list(self) and not ignore_already_added:
                raise EntryAlreadyAddedError(database_entry)
            cmd = commands.AddEntry(self.session, database_entry)
            if self._enable_history:
                cmds.add(cmd)
            else:
                cmd()
            self._cache.append(database_entry)
        if cmds:
            self._command_manager.do(cmds)

    def add_from_file(self, file, ignore_already_added=False):
        """Generate as many database entries as there are FITS headers in the
        given file and add them to the database.

        Parameters
        ----------
        file : str or file-like object
            Either a path pointing to a FITS file or a an opened file-like
            object. If an opened file object, its mode must be one of the
            following rb, rb+, or ab+.

        ignore_already_added : bool, optional
            See :meth:`sunpy.database.Database.add`.

        """
        self.add_many(
            tables.entries_from_file(file, self.default_waveunit),
            ignore_already_added)

    def edit(self, database_entry, **kwargs):
        """Change the given database entry so that it interprets the passed
        key-value pairs as new values where the keys represent the attributes
        of this entry. If no keywords arguments are given, :exc:`ValueError` is
        raised.

        """
        cmd = commands.EditEntry(database_entry, **kwargs)
        if self._enable_history:
            self._command_manager.do(cmd)
        else:
            cmd()
        self._cache[database_entry.id] = database_entry

    def remove_many(self, database_entries):
        """Remove a row of database entries "at once". If this method is used,
        only one entry is saved in the undo history.

        Parameters
        ----------
        database_entries : iterable of sunpy.database.tables.DatabaseEntry
            The database entries that will be removed from the database.
        """
        cmds = CompositeOperation()
        for database_entry in database_entries:
            cmd = commands.RemoveEntry(self.session, database_entry)
            if self._enable_history:
                cmds.add(cmd)
            else:
                cmd()
            try:
                del self._cache[database_entry.id]
            except KeyError:
                pass

        if cmds:
            self._command_manager.do(cmds)

    def remove(self, database_entry):
        """Remove the given database entry from the database table."""
        remove_entry_cmd = commands.RemoveEntry(self.session, database_entry)
        if self._enable_history:
            self._command_manager.do(remove_entry_cmd)
        else:
            remove_entry_cmd()
        try:
            del self._cache[database_entry.id]
        except KeyError:
            # entry cannot be removed because it was already removed or never
            # existed in the database. This can be safely ignored, the user
            # doesn't even know there's a cache here
            pass

    def clear(self):
        """Remove all entries from the database. This operation can be undone
        using the :meth:`undo` method.

        """
        cmds = CompositeOperation()
        for entry in self:
            for tag in entry.tags:
                cmds.add(commands.RemoveTag(self.session, entry, tag))
            # TODO: also remove all FITS header entries and all FITS header
            # comments from each entry before removing the entry itself!
        # remove all entries from all helper tables
        database_tables = [
            tables.JSONDump, tables.Tag, tables.FitsHeaderEntry,
            tables.FitsKeyComment]
        for table in database_tables:
            for entry in self.session.query(table):
                cmds.add(commands.RemoveEntry(self.session, entry))
        for entry in self:
            cmds.add(commands.RemoveEntry(self.session, entry))
            del self._cache[entry.id]
        if self._enable_history:
            self._command_manager.do(cmds)
        else:
            cmds()

    def clear_histories(self):
        """Clears all entries from the undo and redo history.

        See Also
        --------
        :meth:`sunpy.database.commands.CommandManager.clear_histories`
        """
        self._command_manager.clear_histories()  # pragma: no cover

    def undo(self, n=1):
        """undo the last n commands.

        See Also
        --------
        :meth:`sunpy.database.commands.CommandManager.undo`

        """
        self._command_manager.undo(n)  # pragma: no cover

    def redo(self, n=1):
        """redo the last n commands.

        See Also
        --------
        :meth:`sunpy.database.commands.CommandManager.redo`

        """
        self._command_manager.redo(n)  # pragma: no cover

    def __getitem__(self, key):
        if isinstance(key, slice):
            entries = []
            start = 0 if key.start is None else key.start
            stop = len(self) if key.stop is None else key.stop
            step = 1 if key.step is None else key.step
            for i in xrange(start, stop, step):
                try:
                    entry = self[i]
                except IndexError:
                    break
                else:
                    self._cache[entry.id]
                    entries.append(entry)
            return entries
        # support negative indices
        if key < 0 < abs(key) <= len(self):
            key %= len(self)
        for i, entry in enumerate(self):
            if i == key:
                # "touch" the entry in the cache to intentionally cause
                # possible side-effects
                self._cache[entry.id]
                return entry
        raise IndexError

    def __contains__(self, database_entry):
        """Return True if the given database_entry entry is saved in the
        database, False otherwise.

        """
        (ret,), = self.session.query(
            exists().where(tables.DatabaseEntry.id == database_entry.id))
        return ret

    def __iter__(self):
        """iterate over all database entries that have been saved."""
        return iter(self.session.query(tables.DatabaseEntry))

    def __len__(self):
        """Get the number of rows in the table."""
        return self.session.query(tables.DatabaseEntry).count()
