from __future__ import absolute_import

from sqlalchemy import Column, Integer, Float, String, DateTime, Boolean
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class DatabaseEntry(Base):
    __tablename__ = 'data'

    # FIXME: primary key is data provider + file ID + download_time!
    id = Column(Integer, primary_key=True)
    source = Column(String)
    provider = Column(String)
    physobs = Column(String)
    fileid = Column(String)
    observation_time = Column(DateTime)
    instrument = Column(String)
    size = Column(Float)
    mission = Column(String)  # FIXME: does this info really exist?!
    waveunit = Column(String)
    wavemin = Column(Float)
    wavemax = Column(Float)
    path = Column(String)
    download_time = Column(DateTime)
    starred = Column(Boolean)

    @classmethod
    def from_query_result_block(cls, qr_block):
        wave = qr_block.wave
        return cls(
            source=qr_block.source, provider=qr_block.provider,
            physobs=qr_block.pysobs, fileid=qr_block.fileid,
            observation_time=qr_block.time, instrument=qr_block.instrument,
            size=qr_block.size, waveunit=wave.waveunit, wavemin=wave.wavemin,
            wavemax=wave.wavemax)

    def __repr__(self):
        return '<%s(id %d, data provider %s, fileid %s)>' % (
            self.__class__.__name__, self.id, self.provider, self.fileid)
