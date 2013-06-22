from __future__ import absolute_import

from sqlalchemy import create_engine, exists
from sqlalchemy.orm import sessionmaker

from sunpy.database import commands, tables
from sunpy.database.caching import LRUCache


class EntryAlreadyAddedError(Exception):
    """This exception is raised if a database entry is attempted to be added to
    the database although it was already saved in it.

    """
    def __init__(self, database_entry):
        self.database_entry = database_entry

    def __str__(self):
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

    def __str__(self):
        return (
            'the entry %r is already marked '
            'as starred' % self.database_entry)


class EntryAlreadyUnstarredError(Exception):
    """This exception is raised if the star mark from a database entry is
    attempted to be removed although the entry is not starred.

    """
    def __init__(self, database_entry):
        self.database_entry = database_entry

    def __str__(self):
        return (
            'the entry %r is already not marked '
            'as starred' % self.database_entry)


class Database(object):
    """
    .. seealso::

        sqlalchemy.create_engine
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
        """Initialise the database by creating all necessary tables."""
        metadata = tables.Base.metadata
        metadata.create_all(self._engine, checkfirst=checkfirst)

    def commit(self):
        """Flush pending changes and commit the current transaction. This is a
        shortcut for :meth:`session.commit()`.

        .. seealso::

            sqlalchemy.orm.session.Session.commit for more information and
            possible exceptions that may be thrown
        """
        self.session.commit()

    def get_entry_by_id(self, entry_id):
        """Get a database entry by its unique ID number."""
        return self._cache[entry_id]

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

    def get_starred(self):
        """Get all starred database entries as a generator."""
        return (entry for entry in self if entry.starred)

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

        .. seealso::

            :ref:`sunpy.database.CommandManager.undo`

        """
        self._command_manager.undo(n)  # pragma: no cover

    def redo(self, n=1):
        """redo the last n commands.

        .. seealso::

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
        return iter(self.session.query(tables.DatabaseEntry).all())

    def __len__(self):
        """Get the number of rows in the table."""
        return self.session.query(tables.DatabaseEntry).count()
