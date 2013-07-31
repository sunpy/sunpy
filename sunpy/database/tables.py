from __future__ import absolute_import

from time import strptime, mktime
from datetime import datetime
import fnmatch
import os
from itertools import imap

from sqlalchemy import Column, Integer, Float, String, DateTime, Boolean,\
    Table, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

from sunpy.time import parse_time
from sunpy.io import fits
from sunpy.util import print_table

__all__ = [
    'FitsHeaderEntry', 'Tag', 'DatabaseEntry', 'entries_from_query_result',
    'entries_from_path']

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
    """
    DatabaseEntry()

    The class :class:`DatabaseEntry` represents the main table of the database
    and each instance represents one record that *can* be saved in the
    database.

    Attributes
    ----------
    id : int
        A unique ID number. By default it is None, but automatically set to the
        maximum number plus one when this entry is added to the database.
    source : string
        The source is the name of an observatory or the name of a network of
        observatories.
    provider : string
        The name of the server which provides the retrieved data.
    physobs : string
        A physical observable identifier used by VSO.
    fileid : string
        The file ID is a string defined by the data provider that should point
        to a specific data product. The association of fileid to the specific
        data may change sometime, if the fileid always points to the latest
        calibrated data.
    observation_time_start : datetime
        The date and time when the observation of the data started.
    observation_time_end : datetime
        The date and time when the observation of the data ended.
    instrument : string
        The instrument which was used to observe the data.
    size : float
        The size of the data in kilobytes.
    mission : string
        The name of the mission. Currently, this attribute is not set
        automatically by any method.
    waveunit : string
        The wave unit which has been used to measure the data (e.g. 'Angstrom')
    wavemin : float
        The value of the measured wave length.
    wavemax : float
        This is the same value as ``wavemin``. The value is stored twice,
        because each ``suds.sudsobject.QueryResponseBlock`` which is used by
        the vso package contains both these values.
    path : string
        A local file path where the according FITS file is saved.
    download_time : datetime
        The date and time when the files connected to a query have been
        downloaded. Note: this is not the date and time when this entry has
        been added to a database!
    starred : bool
        Entries can be starred to mark them. By default, this value is False.
    fits_header_entries : list
        A list of ``FitsHeaderEntry`` instances.
    tags : list
        A list of ``Tag`` instances. Use :ref:`sunpy.database.Database.tag` to
        add a new tag or multiple tags to a specific entry.

    """
    __tablename__ = 'data'

    # FIXME: primary key is data provider + file ID + download_time!
    id = Column(Integer, primary_key=True)
    source = Column(String)
    provider = Column(String)
    physobs = Column(String)
    fileid = Column(String, unique=True)
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
    starred = Column(Boolean, default=False)
    fits_header_entries = relationship('FitsHeaderEntry', backref='data')
    tags = relationship('Tag', secondary=association_table, backref='data')

    @classmethod
    def from_query_result_block(cls, qr_block):
        """Make a new :class:`DatabaseEntry` instance from a VSO query result
        block. A query result block is usually not created directly; instead,
        one gets instances of ``suds.sudsobject.QueryResponseBlock`` by
        iterating over a VSO query result.

        Examples
        --------
        >>> from sunpy.net import vso
        >>> from sunpy.database import DatabaseEntry
        >>> client = vso.VSOClient()
        >>> qr = client.query(vso.attrs.Time('2001/1/1', '2001/1/2'), vso.attrs.Instrument('eit'))
        >>> DatabaseEntry.from_query_result_block(qr[0])
        <DatabaseEntry(id None, data provider SDAC, fileid /archive/soho/private/data/processed/eit/lz/2001/01/efz20010101.010014)>

        """
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

    @classmethod
    def from_fits_filepath(cls, path):
        """Make a new :class:`DatabaseEntry` instance by using the method
        :meth:`add_fits_header_entries_from_file`. This classmethod is simply a
        shortcut for the following lines::

            entry = DatabaseEntry()
            entry.add_fits_header_entries_from_file(path)

        """
        entry = cls()
        entry.add_fits_header_entries_from_file(path)
        return entry

    def add_fits_header_entries_from_file(self, fits_filepath):
        """Use the header of a FITS file to add this information to this
        database entry. It will be saved in the attribute
        :attr:`fits_header_entries`. If the key INSTRUME, WAVELNTH or
        DATE-OBS / DATE_OBS is available, the attribute :attr:`instrument`,
        :attr:`wavemin` and :attr:`wavemax` or :attr:`observation_time_start`
        is set, respectively.

        Parameters
        ----------
        fits_filepath : file path or file-like object
            File to get header from.  If an opened file object, its mode
            must be one of the following rb, rb+, or ab+).

        Examples
        --------
        >>> from pprint import pprint
        >>> from sunpy.database import DatabaseEntry
        >>> import sunpy
        >>> entry = DatabaseEntry()
        >>> entry.fits_header_entries
        []
        >>> entry.add_fits_header_entries_from_file(sunpy.RHESSI_EVENT_LIST)
        >>> pprint(entry.fits_header_entries)
        [<FitsHeaderEntry(id None, key 'SIMPLE', value True)>,
         <FitsHeaderEntry(id None, key 'BITPIX', value 8)>,
         <FitsHeaderEntry(id None, key 'NAXIS', value 0)>,
         <FitsHeaderEntry(id None, key 'EXTEND', value True)>,
         <FitsHeaderEntry(id None, key 'DATE', value '2011-09-13T15:37:38')>,
         <FitsHeaderEntry(id None, key 'ORIGIN', value 'RHESSI')>,
         <FitsHeaderEntry(id None, key 'OBSERVER', value 'Unknown')>,
         <FitsHeaderEntry(id None, key 'TELESCOP', value 'RHESSI')>,
         <FitsHeaderEntry(id None, key 'INSTRUME', value 'RHESSI')>,
         <FitsHeaderEntry(id None, key 'OBJECT', value 'Sun')>,
         <FitsHeaderEntry(id None, key 'DATE_OBS', value '2002-02-20T11:06:00.000')>,
         <FitsHeaderEntry(id None, key 'DATE_END', value '2002-02-20T11:06:43.330')>,
         <FitsHeaderEntry(id None, key 'TIME_UNI', value 1)>,
         <FitsHeaderEntry(id None, key 'ENERGY_L', value 25.0)>,
         <FitsHeaderEntry(id None, key 'ENERGY_H', value 40.0)>,
         <FitsHeaderEntry(id None, key 'TIMESYS', value '1979-01-01T00:00:00')>,
         <FitsHeaderEntry(id None, key 'TIMEUNIT', value 'd')>]

        """
        # FIXME: store a list of headers and not only the first one!
        header = fits.get_header(fits_filepath)[0]
        self.fits_header_entries.extend(
            FitsHeaderEntry(key, value) for key, value in header.iteritems())
        for header_entry in self.fits_header_entries:
            key, value = header_entry.key, header_entry.value
            if key == 'INSTRUME':
                self.instrument = value
            elif key == 'WAVELNTH':
                self.wavemin = self.wavemax = float(value)
            # NOTE: the key DATE-END or DATE_END is not part of the official
            # FITS standard, but many FITS files use it in their header
            elif key in ('DATE-END', 'DATE_END'):
                self.observation_time_end = parse_time(value)
            elif key in ('DATE-OBS', 'DATE_OBS'):
                self.observation_time_start = parse_time(value)

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
            bool(self.starred) == bool(other.starred) and
            self.fits_header_entries == other.fits_header_entries and
            self.tags == other.tags)

    def __ne__(self, other):  # pragma: no cover
        return not (self == other)

    def __repr__(self):
        return '<%s(id %s, data provider %s, fileid %s)>' % (
            self.__class__.__name__, self.id, self.provider, self.fileid)


def entries_from_query_result(qr):
    """Use a query response returned from :meth:`sunpy.net.vso.VSOClient.query`
    or :meth:`sunpy.net.vso.VSOClient.query_legacy` to generate instances of
    :class:`DatabaseEntry`. Return an iterator over those instances.

    Examples
    --------
    >>> from sunpy.net import vso
    >>> from sunpy.database import entries_from_query_result
    >>> client = vso.VSOClient()
    >>> qr = client.query(vso.attrs.Time('2001/1/1', '2001/1/2'), vso.attrs.Instrument('eit'))
    >>> entries = entries_from_query_result(qr)
    >>> entries.next()
    <DatabaseEntry(id None, data provider SDAC, fileid /archive/soho/private/data/processed/eit/lz/2001/01/efz20010101.010014)>

    """
    return (DatabaseEntry.from_query_result_block(block) for block in qr)


def entries_from_path(fitsdir, recursive=False, pattern='*.fits'):
    """Search the given directory for FITS files and use the corresponding FITS
    headers to generate instances of :class:`DatabaseEntry`.

    Parameters
    ----------
    fitsdir : string
        The directory where to look for FITS files.

    recursive : bool, optional
        If True, the given directory will be searched recursively. Otherwise,
        only the given directory and no subdirectories are searched. The
        default is `False`, i.e. the given directory is not searched
        recursively.

    pattern : string, optional
        The pattern defines how FITS files are detected. The default is to
        collect all files with the filename extension `.fits`. This value is
        passed to the function :func:`fnmatch.filter`, see its documentation
        for more information on the supported syntax.

    Returns
    -------
    generator of (DatabaseEntry, str) pairs
        A generator where each item is a tuple consisting of a
        :class:`DatabaseEntry` instance and the absolute path to the filename
        which was used to make the database entry.

    Examples
    --------
    >>> from pprint import pprint
    >>> from sunpy.database import entries_from_path
    >>> from sunpy.data.test import rootdir as fitsdir
    >>> entries = list(entries_from_path(fitsdir))
    >>> len(entries)
    2
    >>> # and now search `fitsdir` recursive
    >>> entries = list(entries_from_path(fitsdir, True))
    >>> len(entries)
    15
    >>> # print the first 5 items of the FITS header of the first found file
    >>> first_entry, filename = entries[0]
    >>> pprint(first_entry.fits_header_entries[:5])
    [<FitsHeaderEntry(id None, key 'SIMPLE', value True)>,
     <FitsHeaderEntry(id None, key 'BITPIX', value -64)>,
     <FitsHeaderEntry(id None, key 'NAXIS', value 2)>,
     <FitsHeaderEntry(id None, key 'NAXIS1', value 128)>,
     <FitsHeaderEntry(id None, key 'NAXIS2', value 128)>]

    """
    for dirpath, dirnames, filenames in os.walk(fitsdir):
        filename_paths = (os.path.join(dirpath, name) for name in filenames)
        for path in fnmatch.filter(filename_paths, pattern):
            yield DatabaseEntry.from_fits_filepath(path), path
        if not recursive:
            break


# FIXME: document me!
def display_entries(database_entries):
    header = [[
        'ID', 'Source', 'Provider', 'Physobs', 'File ID',
        'Obs. time (start, end)', 'Instrument',
        'Size', 'Wave unit', 'Wave (min, max)',
        'Path', 'Downloaded', 'Starred', 'Tags']]
    rulers = [['=' * len(col) for col in header[0]]]
    data = []
    # make sure that there is no value of type NoneType, hence the many
    # str(...) calls
    for entry in database_entries:
        obs_start_end = '%s %s' % (
            entry.observation_time_start.strftime('%Y%m%dT%H%M%S') or 'N/A',
            entry.observation_time_end.strftime('%Y%m%dT%H%M%S') or 'N/A')
        data.append([
            str(entry.id or 'N/A'),
            str(entry.source or 'N/A'),
            str(entry.provider or 'N/A'),
            str(entry.physobs or 'N/A'),
            str(entry.fileid or 'N/A'),  # TODO: truncate in a sensible way
            obs_start_end,
            str(entry.instrument or 'N/A'),
            str(entry.size or 'N/A'),
            #str(entry.mission),
            str(entry.waveunit or 'N/A'),
            ', '.join(map(str, (entry.wavemin, entry.wavemax))) or 'N/A',
            str(entry.path or 'N/A'),
            str(entry.download_time or 'N/A'),
            'Yes' if entry.starred else 'No',
            ', '.join(imap(str, entry.tags))])
    if not data:
        raise TypeError('given iterable is empty')
    return print_table(header + rulers + data)
