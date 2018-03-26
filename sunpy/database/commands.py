# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

from __future__ import absolute_import

from abc import ABCMeta, abstractmethod
import os

from sqlalchemy.orm import make_transient
from sqlalchemy.exc import InvalidRequestError

from sunpy.extern import six

from sunpy.extern.six.moves import range

__all__ = [
    'EmptyCommandStackError', 'NoSuchEntryError', 'NonRemovableTagError',
    'DatabaseOperation', 'AddEntry', 'RemoveEntry', 'EditEntry',
    'CommandManager']


class EmptyCommandStackError(Exception):
    """This exception is raised if it is attempted to pop from a command stack
    even though it is empty.

    """


class NoSuchEntryError(Exception):
    """This exception is raised if it is attempted to remove an entry even
    though it does not exist in the database.

    """
    def __init__(self, database_entry):
        self.database_entry = database_entry

    def __str__(self):  # pragma: no cover
        return (
            'the database entry {0!r} cannot be removed because it '
            'is not stored in the database'.format(self.database_entry))


class NonRemovableTagError(Exception):
    """This exception is raised if it is attempted to remove a tag from a
    database entry even though it is not saved in this entry.

    """
    def __init__(self, database_entry, tag):
        self.database_entry = tag
        self.tag = tag

    def __str__(self):  # pragma: no cover
        errmsg = 'the tag {0} cannot be removed from the database entry {1!r}'
        return errmsg.format(self.database_entry, self.tag)


@six.add_metaclass(ABCMeta)
class DatabaseOperation(object):
    """This is the abstract main class for all database operations. To
    implement a new operation, inherit from this class and override the methods
    __call__ and undo. Both these methods get no parameters (except for self of
    course). The undo method is expected to do the exact opposite of the
    __call__ method, so that calling __call__ *and* undo multiple times in a
    row must not have any side-effects. This is not checked in any way, though.

    """

    @abstractmethod
    def __call__(self):
        return  # pragma: no cover

    @abstractmethod
    def undo(self):
        return  # pragma: no cover


class CompositeOperation(DatabaseOperation):
    def __init__(self, operations=None):
        if operations is None:
            self._operations = []
        else:
            self._operations = operations

    @property
    def operations(self):
        return self._operations

    def add(self, operation):
        self._operations.append(operation)

    def remove(self, operation):
        self._operations.remove(operation)

    def __call__(self):
        for operation in self._operations:
            # FIXME: What follows is the worst hack of my life. Enjoy.
            # Without it, the test test_clear_database would fail.
            f = open(os.devnull, 'w'); f.write(repr(operation)); f.flush(); f.close()
            operation()

    def undo(self):
        for operation in self._operations:
            operation.undo()

    def __len__(self):
        return len(self._operations)


class AddEntry(DatabaseOperation):
    """Add a new database entry to this session. It is not checked whether an
    equivalent entry is already saved in the session; this has to be checked by
    the caller. The ``undo`` method removes the entry from the session again.

    """
    def __init__(self, session, database_entry):
        self.session = session
        self.database_entry = database_entry

    def __call__(self):
        try:
            self.session.add(self.database_entry)
        except InvalidRequestError:
            # database entry cannot be added because it was removed from the
            # database -> use make_transient to send this object back to
            # the transient state
            make_transient(self.database_entry)
            self.session.add(self.database_entry)

    def undo(self):
        try:
            self.session.delete(self.database_entry)
        except InvalidRequestError:
            # database entry cannot be removed because the last call was not
            # followed by a commit -> use make_transient to revert putting the
            # entry into the pending state
            make_transient(self.database_entry)

    def __repr__(self):
        return '<{0}(session {1!r}, entry id {2})>'.format(
            self.__class__.__name__, self.session, self.database_entry.id)


class RemoveEntry(DatabaseOperation):
    """Remove the given database entry from the session. If it cannot be
    removed, because it is not stored in the session,
    :exc:`sunpy.database.NoSuchEntryError` is raised. The ``undo`` method puts
    the database entry back into the session object.

    """
    def __init__(self, session, entry):
        self.session = session
        self.entry = entry

    def __call__(self):
        try:
            self.session.delete(self.entry)
        except InvalidRequestError:
            # self.database_entry cannot be removed because it's not stored in
            # the database
            raise NoSuchEntryError(self.entry)

    def undo(self):
        make_transient(self.entry)
        self.session.add(self.entry)

    def __repr__(self):
        return '<{0}(session {1!r}, entry {2!r})>'.format(
            self.__class__.__name__, self.session, self.entry)


class EditEntry(DatabaseOperation):
    """Change the properties of the database entry. The given keyword arguments
    are used to set the attributes of the entry. The keys represent the
    attribute name and the values represent the new value of this attribute.
    Example: ``EditEntry(entry, foo='bar')`` will set the attribute ``foo`` of
    ``entry`` to the value ``'bar'``.

    """
    def __init__(self, database_entry, **kwargs):
        self.database_entry = database_entry
        if not kwargs:
            raise ValueError("at least one keyword argument must be given")
        self.kwargs = kwargs
        self.prev_values = {}

    def __call__(self):
        for k, v in six.iteritems(self.kwargs):
            # save those values in the dict prev_values that will be changed
            # so that they can be recovered
            self.prev_values[k] = getattr(self.database_entry, k)
            setattr(self.database_entry, k, v)

    def undo(self):
        for k, v in six.iteritems(self.prev_values):
            setattr(self.database_entry, k, v)

    def __repr__(self):
        return '<EditEntry(kwargs {0!r}, entry id {1})>'.format(
            self.kwargs, self.database_entry.id)


class AddTag(DatabaseOperation):
    def __init__(self, session, database_entry, tag):
        self.session = session
        self.database_entry = database_entry
        self.tag = tag

    def __call__(self):
        try:
            self.database_entry.tags.append(self.tag)
        except InvalidRequestError:
            # self.tag cannot be added because it was just removed
            # -> put it back to transient state
            make_transient(self.tag)
            self.database_entry.tags.append(self.tag)

    def undo(self):
        self.database_entry.tags.remove(self.tag)
        if not self.tag.data:
            # remove the tag from the database as well if it was the last tag
            # assigned to an entry
            try:
                RemoveEntry(self.session, self.tag)()
            except NoSuchEntryError:
                # entry cannot be removed because tag is only connected to
                # entries which are not saved in the database
                # -> can be safely ignored
                pass

    def __repr__(self):
        return "<AddTag(tag '{0}', session {1!r}, entry id {2})>".format(
            self.tag, self.session, self.database_entry.id)


class RemoveTag(DatabaseOperation):
    """Remove the tag from the given database entry. If the tag cannot be
    removed from the database entry because it is not assigned to the entry,
    :exc:`sunpy.database.NonRemovableTagError` is raised. The ``undo`` method
    puts the removed tag back into the tag list of the database entry.

    """
    def __init__(self, session, database_entry, tag):
        self.session = session
        self.database_entry = database_entry
        self.tag = tag

    def __call__(self):
        try:
            self.database_entry.tags.remove(self.tag)
        except ValueError:
            # tag not saved in entry
            raise NonRemovableTagError(self.database_entry, self.tag)
        else:
            if not self.tag.data:
                # remove the tag from the database as well if it was the last tag
                # assigned to an entry
                try:
                    RemoveEntry(self.session, self.tag)()
                except NoSuchEntryError:
                    # entry cannot be removed because tag is only connected to
                    # entries which are not saved in the database
                    # -> can be safely ignored
                    pass

    def undo(self):
        try:
            self.database_entry.tags.append(self.tag)
        except InvalidRequestError:
            # self.tag cannot be added because it was just removed
            # -> put it back to transient state
            try:
                make_transient(self.tag)
                self.database_entry.tags.append(self.tag)
            except InvalidRequestError:
                # self.database_entry has been removed
                # -> put it back to transient state
                make_transient(self.database_entry)
                self.database_entry.tags.append(self.tag)

    def __repr__(self):
        return "<RemoveTag(tag '{0}', session {1!r}, entry id {2})>".format(
            self.tag, self.session, self.database_entry.id)


class CommandManager(object):
    """The CommandManager saves all executed and reverted commands to act as an
    undo-redo-manager. All executed commands are saved in the list attribute
    ``undo_commands`` and all undone commands are saved in the list attribute
    ``redo_commands``. It is not recommended to alter these stacks directly;
    instead, use the methods ``push_undo_command``, ``pop_undo_command``,
    ``push_redo_command``, and ``pop_redo_command``, respectively.

    """
    def __init__(self):
        self.undo_commands = []
        self.redo_commands = []

    def clear_histories(self):
        """Clears all entries from the undo and redo history. If one or
        both of the histories are already empty, no exception is raised.

        """
        del self.undo_commands[:]
        del self.redo_commands[:]

    def push_undo_command(self, command):
        """Push the given command to the undo command stack."""
        self.undo_commands.append(command)

    def pop_undo_command(self):
        """Remove the last command from the undo command stack and return it.
        If the command stack is empty,
        :exc:`sunpy.database.commands.EmptyCommandStackError` is raised.

        """
        try:
            last_undo_command = self.undo_commands.pop()
        except IndexError:
            raise EmptyCommandStackError()
        return last_undo_command

    def push_redo_command(self, command):
        """Push the given command to the redo command stack."""
        self.redo_commands.append(command)

    def pop_redo_command(self):
        """Remove the last command from the redo command stack and return it.
        If the command stack is empty,
        :exc:`sunpy.database.commands.EmptyCommandStackError` is raised.

        """
        try:
            last_redo_command = self.redo_commands.pop()
        except IndexError:
            raise EmptyCommandStackError()
        return last_redo_command

    def do(self, command):
        """Execute the given command (a subclass of DatabaseOperation).
        Exceptions raised from the command are not caught. The passed argument
        may also be an iterable of commands. In this case, every command of the
        iterable is executed and only one entry is saved in the undo history.

        """
        command()
        self.push_undo_command(command)
        # clear the redo stack when a new command was executed
        self.redo_commands[:] = []

    def undo(self, n=1):
        """Undo the last n commands. The default is to undo only the last
        command. If there is no command that can be undone because n is too big
        or because no command has been executed yet,
        :exc:`sunpy.database.commands.EmptyCommandStackError` is raised.

        """
        for _ in range(n):
            command = self.pop_undo_command()
            command.undo()
            self.push_redo_command(command)

    def redo(self, n=1):
        """Redo the last n commands which have been undone using the undo
        method. The default is to redo only the last command which has been
        undone using the undo method. If there is no command that can be redone
        because n is too big or because no command has been undone yet,
        :exc:`sunpy.database.commands.EmptyCommandStackError` is raised.

        """
        for _ in range(n):
            command = self.pop_redo_command()
            command()
            self.push_undo_command(command)
