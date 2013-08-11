# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

from __future__ import absolute_import

from abc import ABCMeta, abstractmethod
import collections

from sqlalchemy.orm import make_transient
from sqlalchemy.exc import InvalidRequestError


__all__ = [
    'EmptyCommandStackError', 'NoSuchEntryError', 'DatabaseOperation',
    'AddEntry', 'RemoveEntry', 'EditEntry', 'CommandManager']


class EmptyCommandStackError(Exception):
    pass


class NoSuchEntryError(Exception):
    def __init__(self, database_entry):
        self.database_entry = database_entry

    def __str__(self):  # pragma: no cover
        return (
            'the database entry %r cannot be removed because it '
            'is not stored in the database' % self.database_entry)


class DatabaseOperation(object):
    """This is the abstract main class for all database operations. To
    implement a new operation, inherit from this class and override the methods
    __call__ and undo. Both these methods get no parameters (except for self of
    course). The undo method is expected to do the exact opposite of the
    __call__ method, so that calling __call__ *and* undo multiple times in a
    row must not have any side-effects. This is not checked in any way, though.

    """
    __metaclass__ = ABCMeta

    def __init__(self, session, database_entry):
        self.session = session
        self.database_entry = database_entry

    @abstractmethod
    def __call__(self):
        return  # pragma: no cover

    @abstractmethod
    def undo(self):
        return  # pragma: no cover


class AddEntry(DatabaseOperation):
    def __call__(self):
        try:
            self.session.add(self.database_entry)
        except InvalidRequestError:
            # database entry cannot be added because it was removed from the
            # database -> use the make_transient to send this object back to
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


class RemoveEntry(DatabaseOperation):
    def __call__(self):
        try:
            self.session.delete(self.database_entry)
        except InvalidRequestError:
            # self.database_entry cannot be removed becaused it's not stored in
            # the database
            raise NoSuchEntryError(self.database_entry)

    def undo(self):
        make_transient(self.database_entry)
        self.session.add(self.database_entry)


class EditEntry(DatabaseOperation):
    def __init__(self, database_entry, **kwargs):
        self.database_entry = database_entry
        if not kwargs:
            raise ValueError("at least one keyword argument must be given")
        self.kwargs = kwargs
        self.prev_values = {}

    def __call__(self):
        for k, v in self.kwargs.iteritems():
            # save those values in the dict prev_values that will be changed
            # so that they can be recovered
            self.prev_values[k] = getattr(self.database_entry, k)
            setattr(self.database_entry, k, v)

    def undo(self):
        for k, v in self.prev_values.iteritems():
            setattr(self.database_entry, k, v)


class RemoveTag(DatabaseOperation):
    def __init__(self, database_entry, tag):
        self.database_entry = database_entry
        self.tag = tag

    def __call__(self):
        self.database_entry.tags.remove(self.tag)

    def undo(self):
        self.database_entry.tags.append(self.tag)


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

    def push_undo_command(self, command):
        """Push the given command to the undo command stack."""
        self.undo_commands.append(command)

    def pop_undo_command(self):
        """Remove the last command from the undo command stack and return it.
        If the command stack is empty, EmptyCommandStackError is raised.

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
        If the command stack is empty, EmptyCommandStackError is raised.

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
        if isinstance(command, collections.Iterable):
            for cmd in command:
                cmd()
        else:
            command()
        self.push_undo_command(command)
        # clear the redo stack when a new command was executed
        self.redo_commands[:] = []

    def undo(self, n=1):
        """Undo the last n commands. The default is to undo only the last
        command. If there is no command that can be undone because n is too big
        or because no command has been executed yet, EmptyCommandStackError is
        raised.

        """
        for _ in xrange(n):
            command = self.pop_undo_command()
            if isinstance(command, collections.Iterable):
                for cmd in command:
                    cmd.undo()
            else:
                command.undo()
            self.push_redo_command(command)

    def redo(self, n=1):
        """Redo the last n commands which have been undone using the undo
        method. The default is to redo only the last command which has been
        undone using the undo method. If there is no command that can be redone
        because n is too big or because no command has been undone yet,
        EmptyCommandStackError is raised.

        """
        for _ in xrange(n):
            command = self.pop_redo_command()
            if isinstance(command, collections.Iterable):
                for cmd in command:
                    cmd()
            else:
                command()
            self.push_undo_command(command)
