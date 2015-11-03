# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).

from __future__ import absolute_import

import fnmatch 

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import pytest

from sunpy.database.commands import AddEntry, RemoveEntry, EditEntry,\
    AddTag, RemoveTag, NoSuchEntryError, NonRemovableTagError,\
    EmptyCommandStackError, CommandManager, CompositeOperation
from sunpy.database.tables import DatabaseEntry, Tag


@pytest.fixture
def session():
    # always create an in-memory database with its own new table in each test
    engine = create_engine('sqlite:///:memory:')
    Session = sessionmaker()
    DatabaseEntry.metadata.create_all(bind=engine)
    return Session(bind=engine)


@pytest.fixture
def command_manager():
    return CommandManager()


def test_add_entry_repr(session):
    entry = DatabaseEntry(id=5)
    repr_result = repr(AddEntry(session, entry))
    expected_repr_result = (
        '<AddEntry('
            'session <sqlalchemy.orm.session.Session object at *>, '
            'entry id 5)>'.format(id(session)))
    assert fnmatch.fnmatch(repr_result, expected_repr_result)


def test_add_entry(session):
    assert not session.new
    entry = DatabaseEntry()
    AddEntry(session, entry)()
    assert len(session.new) == 1
    assert entry.id is None
    session.commit()
    assert not session.new
    assert entry.id == 1


def test_add_entry_undo(session):
    entry = DatabaseEntry()
    cmd = AddEntry(session, entry)
    cmd()
    assert session.query(DatabaseEntry).count() == 1
    assert entry.id == 1
    cmd.undo()
    assert entry in session.deleted
    assert session.query(DatabaseEntry).count() == 0


def test_add_removed_entry(session):
    entry = DatabaseEntry()
    AddEntry(session, entry)()
    session.commit()
    RemoveEntry(session, entry)()
    session.commit()
    AddEntry(session, entry)()
    session.commit()
    assert not session.new
    assert entry.id == 1


def test_add_entry_undo_precommit(session):
    entry = DatabaseEntry()
    cmd = AddEntry(session, entry)
    cmd()
    cmd.undo()
    session.commit()
    assert session.query(DatabaseEntry).count() == 0


def test_edit_entry_repr():
    entry = DatabaseEntry(id=7)
    expected_repr_result = "<EditEntry(kwargs {'foo': 'bar'}, entry id 7)>"
    assert fnmatch.fnmatch(repr(EditEntry(entry, foo='bar')), expected_repr_result)


def test_edit_entry_invalid(session):
    with pytest.raises(ValueError):
        EditEntry(DatabaseEntry())


def test_edit_entry(session):
    entry = DatabaseEntry()
    session.add(entry)
    session.commit()
    assert entry.id == 1
    EditEntry(entry, id=42)()
    assert entry.id == 42


def test_edit_entry_undo(session):
    entry = DatabaseEntry()
    session.add(entry)
    session.commit()
    cmd = EditEntry(entry, id=42)
    cmd()
    session.commit()
    assert entry.id == 42
    cmd.undo()
    session.commit()
    assert entry.id == 1


def test_remove_entry_repr(session):
    entry = DatabaseEntry(id=3)
    expected_repr_result = (
        '<RemoveEntry('
            'session <sqlalchemy.orm.session.Session object at *>, '
            'entry <DatabaseEntry(id 3)>)>'.format(id(session)))
    assert fnmatch.fnmatch(repr(RemoveEntry(session, entry)), expected_repr_result)


def test_remove_existing_entry(session):
    entry = DatabaseEntry()
    session.add(entry)
    assert session.query(DatabaseEntry).count() == 1
    assert entry.id == 1
    RemoveEntry(session, entry)()
    assert entry in session.deleted
    assert session.query(DatabaseEntry).count() == 0


def test_remove_nonexisting_entry(session):
    with pytest.raises(NoSuchEntryError):
        RemoveEntry(session, DatabaseEntry())()


def test_remove_entry_undo(session):
    entry = DatabaseEntry()
    session.add(entry)
    cmd = RemoveEntry(session, entry)
    session.commit()
    cmd()
    assert session.query(DatabaseEntry).count() == 0
    cmd.undo()
    assert session.query(DatabaseEntry).count() == 1


def test_add_tag_repr(session):
    entry = DatabaseEntry(id=12)
    tag = Tag('spam')
    expected_repr_result = (
        "<AddTag("
            "tag 'spam', "
            "session <sqlalchemy.orm.session.Session object at *>, "
            "entry id 12)>".format(id(session)))
    assert fnmatch.fnmatch(repr(AddTag(session, entry, tag)), expected_repr_result)


def test_add_tag(session):
    tag = Tag('tag')
    entry = DatabaseEntry()
    assert entry.tags == []
    cmd = AddTag(session, entry, tag)
    cmd()
    assert tag in entry.tags


def test_add_removed_tag(session):
    entry = DatabaseEntry()
    tag = Tag('tag')
    entry.tags.append(tag)
    session.add(tag)
    session.commit()
    session.delete(tag)
    AddTag(session, entry, tag)()
    assert tag in entry.tags


def test_add_tag_undo_unsaved_entry(session):
    tag = Tag('tag')
    entry = DatabaseEntry()
    cmd = AddTag(session, entry, tag)
    cmd()
    cmd.undo()
    assert entry.tags == []
    cmd()
    assert tag in entry.tags


def test_remove_tag_repr(session):
    entry = DatabaseEntry(id=8)
    tag = Tag('foo')
    expected_repr_result = (
        "<RemoveTag("
            "tag 'foo', "
            "session <sqlalchemy.orm.session.Session object at *>, "
            "entry id 8)>".format(id(session)))
    assert fnmatch.fnmatch(repr(RemoveTag(session, entry, tag)), expected_repr_result)


def test_remove_nonexisting_tag(session):
    cmd = RemoveTag(session, DatabaseEntry(), Tag('tag'))
    with pytest.raises(NonRemovableTagError):
        cmd()


def test_remove_tag_undo(session):
    tag = Tag('tag')
    entry = DatabaseEntry()
    entry.tags.append(tag)
    session.add(entry)
    session.commit()
    assert tag in entry.tags
    cmd = RemoveTag(session, entry, tag)
    cmd()
    assert tag not in entry.tags
    cmd.undo()
    assert tag in entry.tags


def test_cmd_manager_pop_undo_cmd(session, command_manager):
    cmd = AddEntry(session, DatabaseEntry())
    command_manager.do(cmd)
    popped_cmd = command_manager.pop_undo_command()
    assert popped_cmd == cmd


def test_cmd_manager_pop_undo_cmd_empty_stack(command_manager):
    with pytest.raises(EmptyCommandStackError):
        command_manager.pop_undo_command()


def test_cmd_manager_pop_redo_cmd(command_manager):
    with pytest.raises(EmptyCommandStackError):
        command_manager.pop_redo_command()


def test_cmd_manager_pop_redo_cmd_empty_stack(session, command_manager):
    cmd = AddEntry(session, DatabaseEntry())
    command_manager.do(cmd)
    command_manager.undo()
    popped_cmd = command_manager.pop_redo_command()
    assert popped_cmd == cmd


def test_cmd_manager_redo_stack_empty_after_call(session, command_manager):
    command_manager.do(AddEntry(session, DatabaseEntry()))
    command_manager.do(AddEntry(session, DatabaseEntry()))
    assert len(command_manager.undo_commands) == 2
    session.commit()
    command_manager.undo(2)
    assert not command_manager.undo_commands
    assert len(command_manager.redo_commands) == 2
    command_manager.do(AddEntry(session, DatabaseEntry()))
    assert not command_manager.redo_commands


def test_cmd_manager_redo(session, command_manager):
    assert command_manager.undo_commands == []
    assert command_manager.redo_commands == []
    command_manager.do(AddEntry(session, DatabaseEntry()))
    command_manager.do(AddEntry(session, DatabaseEntry()))
    assert len(command_manager.undo_commands) == 2
    assert command_manager.redo_commands == []
    session.commit()
    command_manager.undo(2)
    assert command_manager.undo_commands == []
    assert len(command_manager.redo_commands) == 2
    command_manager.redo(2)
    assert len(command_manager.undo_commands) == 2
    assert command_manager.redo_commands == []


def test_undo_redo_multiple_cmds_at_once(session, command_manager):
    assert command_manager.undo_commands == []
    command_manager.do(CompositeOperation([
        AddEntry(session, DatabaseEntry()),
        AddEntry(session, DatabaseEntry()),
        AddEntry(session, DatabaseEntry())]))
    assert len(command_manager.undo_commands) == 1
    assert session.query(DatabaseEntry).count() == 3
    command_manager.undo()
    assert command_manager.undo_commands == []
    assert session.query(DatabaseEntry).count() == 0
    command_manager.redo()
    assert command_manager.redo_commands == []
    assert session.query(DatabaseEntry).count() == 3
