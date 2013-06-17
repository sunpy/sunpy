from __future__ import absolute_import

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import pytest

from sunpy.database.commands import AddEntry, RemoveEntry, EditEntry,\
    NoSuchEntryError, CommandManager
from sunpy.database.tables import DatabaseEntry


@pytest.fixture
def session():
    # always create an in-memory database with its own new table in each test
    engine = create_engine('sqlite:///:memory:')
    Session = sessionmaker()
    DatabaseEntry.metadata.create_all(bind=engine)
    return Session(bind=engine)


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


def test_redo_stack_empty_after_call(session):
    manager = CommandManager()
    manager.do(AddEntry(session, DatabaseEntry()))
    manager.do(AddEntry(session, DatabaseEntry()))
    assert len(manager.undo_commands) == 2
    session.commit()
    manager.undo(2)
    assert not manager.undo_commands
    assert len(manager.redo_commands) == 2
    manager.do(AddEntry(session, DatabaseEntry()))
    assert not manager.redo_commands
