from __future__ import absolute_import

import operator

from sqlalchemy import create_engine, exists
from sqlalchemy.orm import sessionmaker

from sunpy.database import commands, tables
from sunpy.database.caching import LRUCache
from sunpy.database.attrs import walker
from sunpy.net.attr import and_


class EntryAlreadyAddedError(Exception):
    """This exception is raised if a database entry is attempted to be added to
    the database although it was already saved in it.

    """
    def __init__(self, database_entry):
        self.database_entry = database_entry

    def __str__(self):  # pragma: no cover
        return (
            'the entry %r was already added '
            'to the database' % self.database_entry)


class EntryAlreadyStarredError(Exception):
    """This exception is raised if a database entry is marked as starred
    using :meth:`Database.star` although it was already starred before this
    operation.

    """
    def __init__(self, database_entry):
        self.database_entry = database_entry

    def __str__(self):  # pragma: no cover
        return (
            'the entry %r is already marked '
            'as starred' % self.database_entry)


class EntryAlreadyUnstarredError(Exception):
    """This exception is raised if the star mark from a database entry is
    attempted to be removed although the entry is not starred.

    """
    def __init__(self, database_entry):
        self.database_entry = database_entry

    def __str__(self):  # pragma: no cover
        return (
            'the entry %r is already not marked '
            'as starred' % self.database_entry)


class Database(object):
    """
    Database(url[, CacheClass[, cache_size]])

    Parameters
    ----------
    url : str
        A URL describing the database. This value is simply passed to
        :func:`sqlalchemy.create_engine`
    CacheClass : sunpy.database.caching.BaseCache
        A concrete cache implementation of the abstract class BaseCache.
        Builtin supported values for this parameters are
        :class:`sunpy.database.caching.LRUCache` and
        :class:`sunpy.database.caching.LFUCache`.
        The default value is :class:`sunpy.database.caching.LRUCache`.
    cache_size : int
        The maximum number of database entries, default is no limit.

    Attributes
    ----------
    session : sqlalchemy.orm.session.Session
        A SQLAlchemy session object. This may be used for advanced queries and
        advanced manipulations and should only be used by people who are
        experienced with SQLAlchemy.

    tags : list of sunpy.database.Tag objects
        A list of all saved tags in database. This attribute is read-only.

    Methods
    -------
    create_tables(checkfirst=True)
        Create all necessary tables. Do nothing if ``checkfirst`` is True and
        the required tables already exist.
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
    get_by_tags(*tags)
        Get all database entries that have at least one of the tags
        assigned. If no tags are given, TypeError is raised.
    star(entry, ignore_already_starred=False)
        Mark the given database entry as starred. If ``ignore_already_starred``
        is False and the given entry is already marked as starred,
        EntryAlreadyStarredError is raised.
    unstar(entry, ignore_already_unstarred=False)
        Remove the starred mark of the given entry. If
        ``ignore_already_unstarred`` is False and the entry is not marked as
        starred, EntryAlreadyUnstarredError is raised.
    get_starred()
        Return an iterator over all starred database entries.
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

    Examples
    --------
    >>> database = Database('sqlite:///:memory:', LFUCache, 3)
    >>> database.create_tables()
    >>> database.add(DatabaseEntry())
    >>> database.add(DatabaseEntry())
    >>> database.add(DatabaseEntry())
    >>> pprint(list(database))
    [<DatabaseEntry(id 1, data provider None, fileid None)>,
     <DatabaseEntry(id 2, data provider None, fileid None)>,
     <DatabaseEntry(id 3, data provider None, fileid None)>]
    >>> database.get_entry_by_id(1)
    <DatabaseEntry(id 1, data provider None, fileid None)>
    >>> database.get_entry_by_id(3)
    <DatabaseEntry(id 3, data provider None, fileid None)>
    >>> database.get_entry_by_id(3)
    <DatabaseEntry(id 3, data provider None, fileid None)>
    >>> database.add(DatabaseEntry())
    >>> pprint(list(database))
    [<DatabaseEntry(id 1, data provider None, fileid None)>,
     <DatabaseEntry(id 3, data provider None, fileid None)>,
     <DatabaseEntry(id 4, data provider None, fileid None)>]
    >>> database.add(DatabaseEntry())
    >>> pprint(list(database))
    [<DatabaseEntry(id 1, data provider None, fileid None)>,
     <DatabaseEntry(id 3, data provider None, fileid None)>,
     <DatabaseEntry(id 5, data provider None, fileid None)>]
    >>> database.undo(2)
    >>> pprint(list(database))
    [<DatabaseEntry(id 1, data provider None, fileid None)>,
     <DatabaseEntry(id 3, data provider None, fileid None)>,
     <DatabaseEntry(id 4, data provider None, fileid None)>]

    """
    def __init__(self, url, CacheClass=LRUCache, cache_size=float('inf')):
        self._engine = create_engine(url)
        self._session_cls = sessionmaker(bind=self._engine)
        self.session = self._session_cls()
        self._command_manager = commands.CommandManager()

        class Cache(CacheClass):
            def callback(this, entry_id, database_entry):
                self.remove(database_entry)

            def append(this, value):
                this[max(this or [0]) + 1] = value
        self._cache = Cache(cache_size)

    def create_tables(self, checkfirst=True):
        """Initialise the database by creating all necessary tables. If
        ``checkfirst`` is True, already existing tables are not attempted to be
        created.

        """
        metadata = tables.Base.metadata
        metadata.create_all(self._engine, checkfirst=checkfirst)

    def commit(self):
        """Flush pending changes and commit the current transaction. This is a
        shortcut for :meth:`session.commit()`.

        See Also
        --------
        sqlalchemy.orm.session.Session.commit for more information and
        possible exceptions that may be thrown
        """
        self.session.commit()

    def query(self, *query, **kwargs):
        """
        query(*query[, sortby])
        Send the given query to the database and return a list of
        database entries that satisfy all of the given attributes.

        Apart from the attributes supported by the VSO interface, the following
        attributes are supported:

            - :class:`sunpy.database.attrs.Tag`

            - :class:`sunpy.database.attrs.Starred`

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
            :class:`DatabaseEntry` for a list of all possible values.

        Raises
        ------
        TypeError
            if no attribute is given or if some keyword argument other than
            'sortby' is given.

        Examples
        --------
        The database in this example contains 10 entries, of which the entries
        #4 and #8 have the tag 'foo' and the entries #5 and #10 have the tag
        'bar'; none of them are marked as starred. The query in the following
        example searches for all non-starred entries with the tag 'foo' or
        'bar' (or both).

        >>> from pprint import pprint
        >>> pprint(database.query(~attrs.Starred(), attrs.Tag('foo') | attrs.Tag('bar'), sortby='id'))
        [<DatabaseEntry(id 4, data provider None, fileid None)>,
         <DatabaseEntry(id 5, data provider None, fileid None)>,
         <DatabaseEntry(id 8, data provider None, fileid None)>,
         <DatabaseEntry(id 10, data provider None, fileid None)>]

        """
        if not query:
            raise TypeError('at least one attribute required')
        sortby = kwargs.pop('sortby', 'observation_time_start')
        if kwargs:
            k, v = kwargs.popitem()
            raise TypeError('unexpected keyword argument %r' % k)
        return sorted(
            walker.create(and_(*query), self.session),
            key=operator.attrgetter(sortby))

    def get_entry_by_id(self, entry_id):
        """Get a database entry by its unique ID number."""
        return self._cache[entry_id]

    @property
    def tags(self):
        return self.session.query(tables.Tag).all()

    def get_tag(self, tag_name):
        """Get the tag which has the given name. If no such tag exists, None is
        returned.

        """
        for tag in self.tags:
            if tag_name == tag.name:
                return tag
        return None

    def tag(self, database_entry, *tags):
        """Assign the given database entry the given tags. If no tags are
        given, TypeError is raised.

        """
        if not tags:
            raise TypeError('at least one tag must be given')
        # avoid duplicates
        tag_names = set(tags)
        for tag_name in tag_names:
            tag = self.get_tag(tag_name)
            if tag is None:
                # tag does not exist yet -> create it
                database_entry.tags.append(tables.Tag(tag_name))
            else:
                # tag could be found in the tags table -> add this tag
                database_entry.tags.append(tag)

    def star(self, database_entry, ignore_already_starred=False):
        """Mark the given database entry as starred. If this entry is already
        marked as starred, the behaviour depends on the optional argument
        ``ignore_already_starred``: if it is ``False`` (the default),
        :exc:`sunpy.database.EntryAlreadyStarredError` is raised. Otherwise,
        the entry is kept as starred and no exception is raised.

        """
        if database_entry.starred and not ignore_already_starred:
            raise EntryAlreadyStarredError(database_entry)
        database_entry.starred = True
        self._cache[database_entry.id] = database_entry

    def unstar(self, database_entry, ignore_already_unstarred=False):
        """Remove the starred mark of the given entry. If this entry is not
        marked as starred, the behaviour depends on the optional argument
        ``ignore_already_unstarred``: if it is ``False`` (the default),
        :exc:`sunpy.database.EntryAlreadyUnstarredError` is raised. Otherwise,
        the entry is kept as unstarred and no exception is raised.

        """
        if not database_entry.starred and not ignore_already_unstarred:
            raise EntryAlreadyUnstarredError(database_entry)
        database_entry.starred = False
        self._cache[database_entry.id] = database_entry

    def add(self, database_entry, ignore_already_added=False):
        """Add the given database entry to the database table. If the given
        database entry already exists in the database and
        ``ignore_already_added`` is ``False`` (the default),
        :exc:`sunpy.database.EntryAlreadyAddedError` is raised. Otherwise, the
        entry is kept in the database and no exception is raised.

        """
        if database_entry in self and not ignore_already_added:
            raise EntryAlreadyAddedError(database_entry)
        add_entry_cmd = commands.AddEntry(self.session, database_entry)
        self._command_manager.do(add_entry_cmd)
        if database_entry.id is None:
            self._cache.append(database_entry)
        else:
            self._cache[database_entry.id] = database_entry

    def edit(self, database_entry, **kwargs):
        """Change the given database entry so that it interprets the passed
        key-value pairs as new values where the keys represent the attributes
        of this entry. If no keywords arguments are given, :exc:`ValueError` is
        raised.

        """
        self._command_manager.do(commands.EditEntry(database_entry, **kwargs))
        self._cache[database_entry.id] = database_entry

    def remove(self, database_entry):
        """Remove the given database entry from the database table."""
        remove_entry_cmd = commands.RemoveEntry(self.session, database_entry)
        self._command_manager.do(remove_entry_cmd)
        try:
            del self._cache[database_entry.id]
        except KeyError:
            # entry cannot be removed because it was already removed or never
            # existed in the database. This can be safely ignored, the user
            # doesn't even know there's a cache here
            pass

    def undo(self, n=1):
        """undo the last n commands.

        See Also
        --------
        :ref:`sunpy.database.CommandManager.undo`

        """
        self._command_manager.undo(n)  # pragma: no cover

    def redo(self, n=1):
        """redo the last n commands.

        See Also
        --------
        :ref:`sunpy.database.CommandManager.redo`

        """
        self._command_manager.redo(n)  # pragma: no cover

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
