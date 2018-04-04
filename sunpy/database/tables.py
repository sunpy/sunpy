# Author: Simon Liedtke <liedtke.simon@googlemail.com>
#
# This module was developed with funding provided by
# the Google Summer of Code (2013).
from __future__ import absolute_import, print_function

from time import strptime, mktime
from datetime import datetime
import fnmatch
import os

from astropy.units import Unit, nm, equivalencies, quantity
import astropy.table
from sqlalchemy import Column, Integer, Float, String, DateTime, Boolean,\
    Table, ForeignKey
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
import numpy as np

from sunpy.time import parse_time, TimeRange
from sunpy.io import fits, file_tools as sunpy_filetools
from sunpy.util import print_table
from sunpy.extern.six.moves import map
from sunpy.extern import six
import sunpy.net

from sunpy import config


TIME_FORMAT = config.get("general", "time_format")

__all__ = [
    'WaveunitNotFoundError', 'WaveunitNotConvertibleError', 'JSONDump',
    'FitsHeaderEntry', 'FitsKeyComment', 'Tag', 'DatabaseEntry',
    'entries_from_query_result', 'entries_from_file', 'entries_from_dir',
    'display_entries']

Base = declarative_base()

# required for the many-to-many relation on tags:entries
association_table = Table('association', Base.metadata,
                          Column('tag_name', String, ForeignKey('tags.name')),
                          Column('entry_id', Integer, ForeignKey('data.id'))
                          )


class WaveunitNotFoundError(Exception):
    """This exception is raised if a wavelength unit cannot be found in a FITS
    header or in a VSO query result block.

    """

    def __init__(self, obj):
        self.obj = obj

    def __str__(self):  # pragma: no cover
        return 'the wavelength unit cannot be found in {0}'.format(self.obj) + \
               ' and default_waveunit not specified when opening the database'


class WaveunitNotConvertibleError(Exception):
    """This exception is raised if a wavelength cannot be converted to an
    astropy.units.Unit instance.

    """

    def __init__(self, waveunit):
        self.waveunit = waveunit

    def __str__(self):  # pragma: no cover
        return (
            'the waveunit {0!r} cannot be converted to an '
            'astropy.units.Unit instance'.format(self.waveunit))


# TODO: move this function outside this package (sunpy.util? sunpy.time?)
def timestamp2datetime(format, string):
    return datetime.fromtimestamp(mktime(strptime(string, format)))


class JSONDump(Base):
    __tablename__ = 'jsondumps'

    dump = Column(String, nullable=False, primary_key=True)

    def __init__(self, dump):
        self.dump = dump

    def __eq__(self, other):
        return self.dump == other.dump

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        return self.dump

    def __repr__(self):  # pragma: no cover
        return '<{0}(dump {1!r})>'.format(self.__class__.__name__, self.dump)


class FitsHeaderEntry(Base):
    __tablename__ = 'fitsheaderentries'

    dbentry_id = Column(Integer, ForeignKey('data.id'))
    id = Column(Integer, primary_key=True)
    key = Column(String, nullable=False)
    value = Column(String)

    def __init__(self, key, value):
        self.key = key
        self.value = value

    def __eq__(self, other):
        return (
            (self.id == other.id or self.id is None or other.id is None) and
            self.key == other.key and
            self.value == other.value)

    def __hash__(self):
        return super(FitsHeaderEntry, self).__hash__()

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):  # pragma: no cover
        return '<{0}(id {1}, key {2!r}, value {3!r})>'.format(
            self.__class__.__name__, self.id, self.key, self.value)


class FitsKeyComment(Base):
    __tablename__ = 'fitskeycomments'

    dbentry_id = Column(Integer, ForeignKey('data.id'))
    id = Column(Integer, primary_key=True)
    key = Column(String, nullable=False)
    value = Column(String)

    def __init__(self, key, value):
        self.key = key
        self.value = value

    def __eq__(self, other):
        return (
            (self.id == other.id or self.id is None or other.id is None) and
            self.key == other.key and
            self.value == other.value)

    def __lt__(self, other):
        return ('{0}, {1}'.format(self.key, self.value) <
                '{0}, {1}'.format(other.key, other.value))

    def __hash__(self):
        return super(FitsKeyComment, self).__hash__()

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):  # pragma: no cover
        return '<{0}(id {1}, key {2!r}, value {3!r})>'.format(
            self.__class__.__name__, self.id, self.key, self.value)


class Tag(Base):
    __tablename__ = 'tags'

    name = Column(String, nullable=False, primary_key=True)

    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return super(Tag, self).__hash__()

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        return self.name

    def __repr__(self):  # pragma: no cover
        return '<{0}(name {1!r})>'.format(self.__class__.__name__, self.name)


class DatabaseEntry(Base):
    """
    DatabaseEntry()

    The class :class:`DatabaseEntry` represents the main table of the database
    and each instance represents one record that *can* be saved in the
    database.

    Parameters
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
    wavemin : float
        The value of the measured wave length.
    wavemax : float
        This is the same value as ``wavemin``. The value is stored twice,
        because each ``suds.sudsobject.QueryResponseBlock`` which is used by
        the vso package contains both these values.
    hdu_index : int
        This value provides a list of all available HDUs and in what
        files they are located.
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
        A list of ``Tag`` instances. Use `sunpy.database.Database.tag` to
        add a new tag or multiple tags to a specific entry.

    """
    __tablename__ = 'data'

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
    wavemin = Column(Float)
    wavemax = Column(Float)
    hdu_index = Column(Integer)
    path = Column(String)
    download_time = Column(DateTime)
    starred = Column(Boolean, default=False)
    fits_header_entries = relationship('FitsHeaderEntry')
    fits_key_comments = relationship('FitsKeyComment')
    tags = relationship('Tag', secondary=association_table, backref='data')

    @classmethod
    def _from_query_result_block(cls, qr_block, default_waveunit=None):
        """Make a new :class:`DatabaseEntry` instance from a VSO query result
        block. The values of :attr:`wavemin` and :attr:`wavemax` are converted
        to nm (nanometres).

        Parameters
        ----------
        qr_block : suds.sudsobject.QueryResponseBlock
            A query result block is usually not created directly; instead,
            one gets instances of ``suds.sudsobject.QueryResponseBlock`` by
            iterating over a VSO query result.
        default_waveunit : str, optional
            The wavelength unit that is used if it cannot be found in the
            `qr_block`.

        Examples
        --------
        >>> from sunpy.net import vso
        >>> from sunpy.database.tables import DatabaseEntry
        >>> client = vso.VSOClient()  # doctest: +REMOTE_DATA
        >>> qr = client.search(
        ...     vso.attrs.Time('2001/1/1', '2001/1/2'),
        ...     vso.attrs.Instrument('eit'))  # doctest: +REMOTE_DATA
        >>> entry = DatabaseEntry._from_query_result_block(qr[0])  # doctest: +REMOTE_DATA
        >>> entry.source  # doctest: +REMOTE_DATA
        SOHO
        >>> entry.provider  # doctest: +REMOTE_DATA
        SDAC
        >>> entry.physobs  # doctest: +REMOTE_DATA
        'intensity'
        >>> entry.fileid  # doctest: +REMOTE_DATA
        /archive/soho/private/data/processed/eit/lz/2001/01/efz20010101.000042
        >>> entry.observation_time_start, entry.observation_time_end  # doctest: +REMOTE_DATA
        (datetime.datetime(2001, 1, 1, 0, 0, 42), datetime.datetime(2001, 1, 1, 0, 0, 54))
        >>> entry.instrument  # doctest: +REMOTE_DATA
        EIT
        >>> entry.size  # doctest: +REMOTE_DATA
        2059.0
        >>> entry.wavemin, entry.wavemax  # doctest: +REMOTE_DATA
        (19.5, 19.5)

        """
        time_start = timestamp2datetime('%Y%m%d%H%M%S', qr_block.time.start)
        if not qr_block.time.end:
            qr_block.time.end = qr_block.time.start
        time_end = timestamp2datetime('%Y%m%d%H%M%S', qr_block.time.end)
        wave = qr_block.wave
        unit = None
        if wave.waveunit is None:
            if default_waveunit is not None:
                unit = Unit(default_waveunit)
        else:
            # some query response blocks store the unit "kev",
            # but Astropy only understands "keV". See issue #766.
            waveunit = wave.waveunit
            if waveunit == "kev":
                waveunit = "keV"
            unit = Unit(waveunit)
        if wave.wavemin is None:
            wavemin = None
        else:
            if unit is None:
                raise WaveunitNotFoundError(qr_block)
            wavemin = unit.to(nm, float(wave.wavemin),
                              equivalencies.spectral())
        if wave.wavemax is None:
            wavemax = None
        else:
            if unit is None:
                raise WaveunitNotFoundError(qr_block)
            wavemax = unit.to(nm, float(wave.wavemax),
                              equivalencies.spectral())
        source = getattr(qr_block, 'source', None)
        provider = getattr(qr_block, 'provider', None)
        fileid = getattr(qr_block, 'fileid', None)
        instrument = getattr(qr_block, 'instrument', None)
        size = getattr(qr_block, 'size', -1)
        physobs = getattr(qr_block, 'physobs', None)
        if physobs is not None:
            physobs = str(physobs)
        return cls(
            source=source, provider=provider, physobs=physobs, fileid=fileid,
            observation_time_start=time_start, observation_time_end=time_end,
            instrument=instrument, size=size,
            wavemin=wavemin, wavemax=wavemax)

    @classmethod
    def _from_fido_search_result_block(cls, sr_block, default_waveunit=None):
        """
        Make a new :class:`DatabaseEntry` instance from a Fido search
        result block.

        Parameters
        ----------
        sr_block : `sunpy.net.dataretriever.client.QueryResponseBlock`
            A query result block is usually not created directly; instead,
            one gets instances of
            ``sunpy.net.dataretriever.client.QueryResponseBlock`` by iterating
            over each element of a Fido search result.
        default_waveunit : `str`, optional
            The wavelength unit that is used if it cannot be found in the
            `sr_block`.
        """
        # All attributes of DatabaseEntry that are not in QueryResponseBlock
        # are set as None for now.
        source = getattr(sr_block, 'source', None)
        provider = getattr(sr_block, 'provider', None)
        physobs = getattr(sr_block, 'physobs', None)
        if physobs is not None:
            physobs = str(physobs)
        instrument = getattr(sr_block, 'instrument', None)
        time_start = sr_block.time.start
        time_end = sr_block.time.end

        wavelengths = getattr(sr_block, 'wave', None)
        wavelength_temp = {}
        if isinstance(wavelength_temp, tuple):
            # Tuple of values
            wavelength_temp['wavemin'] = wavelengths[0]
            wavelength_temp['wavemax'] = wavelengths[1]
        else:
            # Single Value
            wavelength_temp['wavemin'] = wavelength_temp['wavemax'] = wavelengths

        final_values = {}
        for key, val in wavelength_temp.items():
            if isinstance(val, quantity.Quantity):
                unit = getattr(val, 'unit', None)
                if unit is None:
                    if default_waveunit is not None:
                        unit = Unit(default_waveunit)
                    else:
                        raise WaveunitNotFoundError(sr_block)
                final_values[key] = unit.to(nm, float(val.value), equivalencies.spectral())
            elif val is None or np.isnan(val):
                final_values[key] = val

        wavemin = final_values['wavemin']
        wavemax = final_values['wavemax']

        # sr_block.url of a QueryResponseBlock attribute is stored in fileid
        fileid = str(sr_block.url) if sr_block.url is not None else None
        size = None
        return cls(
            source=source, provider=provider, physobs=physobs, fileid=fileid,
            observation_time_start=time_start, observation_time_end=time_end,
            instrument=instrument, size=size,
            wavemin=wavemin, wavemax=wavemax)

    def __eq__(self, other):

        if self.wavemin is None and other.wavemin is None:
            wavemins_equal = True
        elif not all([self.wavemin, other.wavemin]):
            # This means one is None and the other isnt
            wavemins_equal = False
        else:
            wavemins_equal = np.allclose([self.wavemin], [other.wavemin], equal_nan=True)

        if self.wavemax is None and other.wavemax is None:
            wavemaxs_equal = True
        elif not all([self.wavemax, other.wavemax]):
            # This means one is None and the other isnt
            wavemaxs_equal = False
        else:
            wavemaxs_equal = np.allclose([self.wavemax], [other.wavemax], equal_nan=True)

        return (
            (self.id == other.id or self.id is None or other.id is None) and
            self.source == other.source and
            self.provider == other.provider and
            self.physobs == other.physobs and
            self.fileid == other.fileid and
            self.observation_time_start == other.observation_time_start and
            self.observation_time_end == other.observation_time_end and
            self.instrument == other.instrument and
            self.size == other.size and
            wavemins_equal and
            wavemaxs_equal and
            self.path == other.path and
            self.download_time == other.download_time and
            bool(self.starred) == bool(other.starred) and
            self.fits_header_entries == other.fits_header_entries and
            self.tags == other.tags)

    def _compare_attributes(self, other, attribute_list):
        """
        Compare a given list of attributes of two :class:`DatabaseEntry`
        instances and return True if all of them match.

        Parameters
        ----------
        other : :class:`DatabaseEntry` instance

        attribute_list : `list`
            The list of attributes that will be compared in both instances,
            self and other.

        """
        if len(attribute_list) == 0:
            raise TypeError('At least one attribute required')
        for attribute in attribute_list:
            if getattr(self, attribute) != getattr(other, attribute):
                return False
        return True

    def __hash__(self):
        return super(DatabaseEntry, self).__hash__()

    def __ne__(self, other):  # pragma: no cover
        return not (self == other)

    def __repr__(self):  # pragma: no cover
        attrs = [
            'id', 'source', 'provider', 'physobs', 'fileid',
            'observation_time_start', 'observation_time_end', 'instrument',
            'size', 'wavemin', 'wavemax', 'path', 'download_time', 'starred',
            'fits_header_entries', 'tags']
        ret = '<{0}('.format(self.__class__.__name__)
        for attr in attrs:
            value = getattr(self, attr, None)
            if value:
                ret += '{0} {1!r}, '.format(attr, value)
        ret = ret.rstrip(', ')
        ret += ')>'
        return ret


def entries_from_query_result(qr, default_waveunit=None):
    """
    Use a query response returned from :meth:`sunpy.net.vso.VSOClient.search`
    or :meth:`sunpy.net.vso.VSOClient.query_legacy` to generate instances of
    :class:`DatabaseEntry`. Return an iterator over those instances.

    Parameters
    ----------
    qr : `sunpy.net.vso.QueryResponse`
        The query response from which to build the database entries.

    default_waveunit : `str`, optional
        See :meth:`sunpy.database.DatabaseEntry.from_query_result_block`.

    Examples
    --------
    >>> from sunpy.net import vso
    >>> from sunpy.database.tables import entries_from_query_result
    >>> client = vso.VSOClient()  # doctest: +REMOTE_DATA
    >>> qr = client.search(
    ...     vso.attrs.Time('2001/1/1', '2001/1/2'),
    ...     vso.attrs.Instrument('eit'))  # doctest: +REMOTE_DATA
    >>> entries = entries_from_query_result(qr)  # doctest: +REMOTE_DATA
    >>> entry = next(entries)  # doctest: +REMOTE_DATA
    >>> entry.source  # doctest: +REMOTE_DATA
    SOHO
    >>> entry.provider  # doctest: +REMOTE_DATA
    SDAC
    >>> entry.physobs  # doctest: +REMOTE_DATA
    'intensity'
    >>> entry.fileid  # doctest: +REMOTE_DATA
    /archive/soho/private/data/processed/eit/lz/2001/01/efz20010101.000042
    >>> entry.observation_time_start, entry.observation_time_end  # doctest: +REMOTE_DATA
    (datetime.datetime(2001, 1, 1, 0, 0, 42), datetime.datetime(2001, 1, 1, 0, 0, 54))
    >>> entry.instrument  # doctest: +REMOTE_DATA
    EIT
    >>> entry.size  # doctest: +REMOTE_DATA
    2059.0
    >>> entry.wavemin, entry.wavemax  # doctest: +REMOTE_DATA
    (19.5, 19.5)

    """
    for block in qr:
        yield DatabaseEntry._from_query_result_block(block, default_waveunit)


def entries_from_fido_search_result(sr, default_waveunit=None):
    """
    Use a `sunpy.net.dataretriever.fido_factory.UnifiedResponse`
    object returned from
    :meth:`sunpy.net.dataretriever.fido_factory.UnifiedDownloaderFactory.search`
    to generate instances of :class:`DatabaseEntry`. Return an iterator
    over those instances.

    Parameters
    ----------
    search_result : `sunpy.net.dataretriever.fido_factory.UnifiedResponse`
            A UnifiedResponse object that is used to store responses from the
            unified downloader. This is returned by the ``search`` method of a
            :class:`sunpy.net.dataretriever.fido_factory.UnifiedDownloaderFactory`
            object.

    default_waveunit : `str`, optional
        The wavelength unit that is used if it cannot be found in the Query
        Response block.

    Examples
    --------
    >>> from sunpy.net import Fido, attrs
    >>> from sunpy.database.tables import entries_from_fido_search_result
    >>> sr = Fido.search(attrs.Time("2012/1/1", "2012/1/2"),
    ...     attrs.Instrument('lyra')) # doctest: +REMOTE_DATA
    >>> entries = entries_from_fido_search_result(sr) # doctest: +REMOTE_DATA
    >>> entry = next(entries) # doctest: +REMOTE_DATA
    >>> entry.source # doctest: +REMOTE_DATA
    'Proba2'
    >>> entry.provider # doctest: +REMOTE_DATA
    'esa'
    >>> entry.physobs # doctest: +REMOTE_DATA
    'irradiance'
    >>> entry.fileid # doctest: +REMOTE_DATA
    'http://proba2.oma.be/lyra/data/bsd/2012/01/01/lyra_20120101-000000_lev2_std.fits'
    >>> entry.observation_time_start, entry.observation_time_end # doctest: +REMOTE_DATA
    (datetime.datetime(2012, 1, 1, 0, 0), datetime.datetime(2012, 1, 2, 0, 0))
    >>> entry.instrument # doctest: +REMOTE_DATA
    'lyra'

    """
    for entry in sr:
        if isinstance(entry, sunpy.net.vso.vso.QueryResponse):
            # This is because Fido can search the VSO. It
            # returns a VSO QueryResponse.
            for block in entry:
                yield DatabaseEntry._from_query_result_block(block, default_waveunit)
        elif isinstance(entry, sunpy.net.jsoc.jsoc.JSOCResponse):
            # Adding JSOC results to the DB not supported for now
            raise ValueError("Cannot add JSOC results to database")
        else:
            for block in entry:
                yield DatabaseEntry._from_fido_search_result_block(block, default_waveunit)


def entries_from_file(file, default_waveunit=None,
                      time_string_parse_format=None):
    """Use the headers of a FITS file to generate an iterator of
    :class:`sunpy.database.tables.DatabaseEntry` instances. Gathered
    information will be saved in the attribute `fits_header_entries`. If the
    key INSTRUME, WAVELNTH or DATE-OBS / DATE_OBS is available, the attribute
    `instrument`, `wavemin` and `wavemax` or `observation_time_start` is set,
    respectively. If the wavelength unit can be read, the values of `wavemin`
    and `wavemax` are converted to nm (nanometres). The value of the `file`
    parameter is used to set the attribute `path` of each generated database
    entry.

    Parameters
    ----------
    file : str or file-like object
        Either a path pointing to a FITS file or a an opened file-like object.
        If an opened file object, its mode must be one of the following rb,
        rb+, or ab+.

    default_waveunit : str, optional
        The wavelength unit that is used for a header if it cannot be
        found.

    time_string_parse_format : str, optional
        Fallback timestamp format which will be passed to
        `~datetime.datetime.strftime` if `sunpy.time.parse_time` is unable to
        automatically read the `date-obs` metadata.

    Raises
    ------
    sunpy.database.WaveunitNotFoundError
        If `default_waveunit` is not given and the wavelength unit cannot
        be found in one of the FITS headers

    sunpy.WaveunitNotConvertibleError
        If a wavelength unit could be found but cannot be used to create an
        instance of the type ``astropy.units.Unit``. This can be the case
        for example if a FITS header has the key `WAVEUNIT` with the value
        `nonsense`.

    Examples
    --------
    >>> from sunpy.database.tables import entries_from_file
    >>> import sunpy.data.sample  # doctest: +REMOTE_DATA
    >>> entries = list(entries_from_file(sunpy.data.sample.SWAP_LEVEL1_IMAGE))  # doctest: +REMOTE_DATA
    >>> len(entries)  # doctest: +REMOTE_DATA
    1
    >>> entry = entries.pop()  # doctest: +REMOTE_DATA
    >>> entry.instrument  # doctest: +REMOTE_DATA
    'SWAP'
    >>> entry.observation_time_start, entry.observation_time_end  # doctest: +REMOTE_DATA
    (datetime.datetime(2011, 6, 7, 6, 33, 29, 759000), None)
    >>> entry.wavemin, entry.wavemax  # doctest: +REMOTE_DATA
    (17.400000000000002, 17.400000000000002)
    >>> len(entry.fits_header_entries)  # doctest: +REMOTE_DATA
    111

    """
    headers = fits.get_header(file)
    if isinstance(file, (str, six.text_type)):
        filename = file
    else:
        filename = getattr(file, 'name', None)
    for header in headers:
        entry = DatabaseEntry(path=filename)
        for key, value in six.iteritems(header):
            # Yes, it is possible to have an empty key in a FITS file.
            # Example: sunpy.data.sample.EIT_195_IMAGE
            # Don't ask me why this could be a good idea.
            if key == '':
                value = str(value)
            elif key == 'KEYCOMMENTS':
                for k, v in six.iteritems(value):
                    entry.fits_key_comments.append(FitsKeyComment(k, v))
                continue
            entry.fits_header_entries.append(FitsHeaderEntry(key, value))
        waveunit = fits.extract_waveunit(header)
        entry.hdu_index = headers.index(header)
        if waveunit is None:
            waveunit = default_waveunit
        unit = None
        if waveunit is not None:
            try:
                unit = Unit(waveunit)
            except ValueError:
                raise WaveunitNotConvertibleError(waveunit)
        for header_entry in entry.fits_header_entries:
            key, value = header_entry.key, header_entry.value
            if key == 'INSTRUME':
                entry.instrument = value
            elif key == 'WAVELNTH':
                if unit is None:
                    raise WaveunitNotFoundError(file)
                # use the value of `unit` to convert the wavelength to nm
                entry.wavemin = entry.wavemax = unit.to(
                    nm, value, equivalencies.spectral())
            # NOTE: the key DATE-END or DATE_END is not part of the official
            # FITS standard, but many FITS files use it in their header
            elif key in ('DATE-END', 'DATE_END'):
                entry.observation_time_end = parse_time(
                    value,
                    _time_string_parse_format=time_string_parse_format
                )
            elif key in ('DATE-OBS', 'DATE_OBS'):
                entry.observation_time_start = parse_time(
                    value,
                    _time_string_parse_format=time_string_parse_format
                )
        yield entry


def entries_from_dir(fitsdir, recursive=False, pattern='*',
                     default_waveunit=None, time_string_parse_format=None):
    """Search the given directory for FITS files and use the corresponding FITS
    headers to generate instances of :class:`DatabaseEntry`. FITS files are
    detected by reading the content of each file, the `pattern` argument may be
    used to avoid reading entire directories if one knows that all FITS files
    have the same filename extension.

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
        The pattern can be used to filter the list of filenames before the
        files are attempted to be read. The default is to collect all files.
        This value is passed to the function :func:`fnmatch.filter`, see its
        documentation for more information on the supported syntax.

    default_waveunit : str, optional
        See
        :meth:`sunpy.database.tables.DatabaseEntry.add_fits_header_entries_from_file`.

    time_string_parse_format : str, optional
        Fallback timestamp format which will be passed to
        `~datetime.datetime.strftime` if `sunpy.time.parse_time` is unable to
        automatically read the `date-obs` metadata.

    Returns
    -------
    generator of (DatabaseEntry, str) pairs
        A generator where each item is a tuple consisting of a
        :class:`DatabaseEntry` instance and the absolute path to the filename
        which was used to make the database entry.

    Examples
    --------
    >>> import os
    >>> from sunpy.data.test import rootdir as fitsdir
    >>> from sunpy.database.tables import entries_from_dir

    >>> eitdir = os.path.join(fitsdir, 'EIT')
    >>> entries = list(entries_from_dir(eitdir, default_waveunit='angstrom'))
    >>> len(entries)
    13

    """
    for dirpath, dirnames, filenames in os.walk(fitsdir):
        filename_paths = (os.path.join(dirpath, name) for name in filenames)
        for path in fnmatch.filter(filename_paths, pattern):
            try:
                filetype = sunpy_filetools._detect_filetype(path)
            except (
                    sunpy_filetools.UnrecognizedFileTypeError,
                    sunpy_filetools.InvalidJPEG2000FileExtension):
                continue
            if filetype == 'fits':
                for entry in entries_from_file(
                        path, default_waveunit,
                        time_string_parse_format=time_string_parse_format
                ):
                    yield entry, path
        if not recursive:
            break


def _create_display_table(database_entries, columns=None, sort=False):
    """Generate a table to display the database entries.

    Parameters
    ----------
    database_entries : iterable of :class:`DatabaseEntry` instances
        The database entries will be the rows in the resulting table.

    columns : iterable of str
        The columns that will be displayed in the resulting table. Possible
        values for the strings are all attributes of :class:`DatabaseEntry`.

    sort : bool (optional)
        If True, sorts the entries before displaying them.

    Returns
    -------
    str
        An astropy table that can be printed on the console or written to a
        file.

    """
    if columns is None:
        columns = ['id', 'observation_time_start', 'observation_time_end',
                   'instrument', 'source', 'provider', 'physobs', 'wavemin',
                   'wavemax', 'path', 'fileid', 'tags', 'starred',
                   'download_time', 'size']

    data = []
    for entry in database_entries:
        row = []
        for col in columns:
            if col == 'starred':
                row.append('Yes' if entry.starred else 'No')
            elif col == 'tags':
                row.append(', '.join(map(str, entry.tags)) or 'N/A')
            elif col == 'hdu_index':
                row.append(entry.hdu_index)
            # do not display microseconds in datetime columns
            elif col in (
                    'observation_time_start',
                    'observation_time_end',
                    'download_time'):
                time = getattr(entry, col, None)
                if time is None:
                    formatted_time = 'N/A'
                else:
                    formatted_time = time.strftime(TIME_FORMAT)
                row.append(formatted_time)
            else:
                row.append(str(getattr(entry, col) or 'N/A'))
        if not row:
            raise TypeError('at least one column must be given')
        data.append(row)
    if not data:
        raise TypeError('given iterable is empty')
    if sort:
        data.sort()
    return astropy.table.Table(rows=data, names=columns)


def display_entries(database_entries, columns=None, sort=False):
    """Print a table to display the database entries.

    Parameters
    ----------
    database_entries : iterable of :class:`DatabaseEntry` instances
        The database entries will be the rows in the resulting table.

    columns : iterable of str
        The columns that will be displayed in the resulting table. Possible
        values for the strings are all attributes of :class:`DatabaseEntry`.

    sort : bool (optional)
        If True, sorts the entries before displaying them.

    """
    return _create_display_table(database_entries, columns, sort).__str__()
