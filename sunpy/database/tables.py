from __future__ import absolute_import

from time import strptime, mktime
from datetime import datetime

from sqlalchemy import Column, Integer, Float, String, DateTime, Boolean,\
    Table, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from pyfits import getheader as get_pyfits_header

__all__ = [
    'FitsHeaderEntry', 'Tag', 'DatabaseEntry', 'entries_from_query_result']

Base = declarative_base()

# required for the many-to-many relation on tags:entries
association_table = Table('association', Base.metadata,
    Column('tag_id', Integer, ForeignKey('tags.id')),
    Column('entry_id', Integer, ForeignKey('data.id'))
)


# TODO: move this function outside this package (sunpy.util? sunpy.time?)
def timestamp2datetime(format, string):
    return datetime.fromtimestamp(mktime(strptime(string, format)))


class FitsHeaderEntry(Base):
    __tablename__ = 'fitsheaderentries'
    __hash__ = None

    id = Column(Integer, primary_key=True)
    dbentry_id = Column(Integer, ForeignKey('data.id'))
    key = Column(String, nullable=False)
    value = Column(String, nullable=False)

    def __init__(self, key, value):
        self.key = key
        self.value = value

    def __eq__(self, other):
        return (
            self.id == other.id and
            self.key == other.key and
            self.value == other.value)

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):  # pragma: no cover
        return '<%s(id %s, key %r, value %r)>' % (
            self.__class__.__name__, self.id, self.key, self.value)


class Tag(Base):
    __tablename__ = 'tags'

    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False, unique=True)

    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        return self.id == other.id and self.name == other.name

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):  # pragma: no cover
        return '<%s(id %s, name %r)>' % (
            self.__class__.__name__, self.id, self.name)


class DatabaseEntry(Base):
    __tablename__ = 'data'
    __hash__ = None

    # FIXME: primary key is data provider + file ID + download_time!
    id = Column(Integer, primary_key=True)
    source = Column(String)
    provider = Column(String)
    physobs = Column(String)
    fileid = Column(String)
    observation_time_start = Column(DateTime)
    observation_time_end = Column(DateTime)
    instrument = Column(String)
    size = Column(Float)
    mission = Column(String)  # FIXME: does this info really exist?!
    waveunit = Column(String)
    wavemin = Column(Float)
    wavemax = Column(Float)
    path = Column(String)
    download_time = Column(DateTime)
    starred = Column(Boolean)
    fits_header_entries = relationship('FitsHeaderEntry', backref='data')
    tags = relationship('Tag', secondary=association_table, backref='data')

    @classmethod
    def from_query_result_block(cls, qr_block):
        time_start = timestamp2datetime('%Y%m%d%H%M%S', qr_block.time.start)
        time_end = timestamp2datetime('%Y%m%d%H%M%S', qr_block.time.end)
        wave = qr_block.wave
        return cls(
            source=qr_block.source, provider=qr_block.provider,
            physobs=qr_block.physobs, fileid=qr_block.fileid,
            observation_time_start=time_start, observation_time_end=time_end,
            instrument=qr_block.instrument, size=qr_block.size,
            waveunit=wave.waveunit, wavemin=float(wave.wavemin),
            wavemax=float(wave.wavemax))

    def add_fits_header_entries_from_file(self, fits_filepath):
        header = get_pyfits_header(fits_filepath)
        fits_header_entries = [
            FitsHeaderEntry(key, value) for key, value in header.iteritems()]
        self.fits_header_entries.extend(fits_header_entries)

    def __eq__(self, other):
        return (
            self.id == other.id and
            self.source == other.source and
            self.provider == other.provider and
            self.physobs == other.physobs and
            self.fileid == other.fileid and
            self.observation_time_start == other.observation_time_start and
            self.observation_time_end == other.observation_time_end and
            self.instrument == other.instrument and
            self.size == other.size and
            self.mission == other.mission and
            self.waveunit == other.waveunit and
            self.wavemin == other.wavemin and
            self.wavemax == other.wavemax and
            self.path == other.path and
            self.download_time == other.download_time and
            self.starred == other.starred and
            self.fits_header_entries == other.fits_header_entries and
            self.tags == other.tags)

    def __ne__(self, other):  # pragma: no cover
        return not (self == other)

    def __repr__(self):
        return '<%s(id %s, data provider %s, fileid %s)>' % (
            self.__class__.__name__, self.id, self.provider, self.fileid)


def entries_from_query_result(qr):
    return (DatabaseEntry.from_query_result_block(block) for block in qr)
