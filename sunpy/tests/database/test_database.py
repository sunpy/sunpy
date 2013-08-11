from __future__ import absolute_import

from datetime import datetime

import pytest
import sqlalchemy

import sunpy
from sunpy.database import Database, EntryAlreadyAddedError,\
    EntryAlreadyStarredError, EntryAlreadyUnstarredError, NoSuchTagError,\
    EntryNotFoundError, TagAlreadyAssignedError
from sunpy.database.tables import DatabaseEntry, FitsHeaderEntry, Tag, display_entries
from sunpy.database.commands import NoSuchEntryError
from sunpy.database.caching import LRUCache, LFUCache
from sunpy.database import attrs
from sunpy.net import vso


@pytest.fixture
def database_without_tables():
    return Database('sqlite:///:memory:')


@pytest.fixture
def database_using_lrucache():
    d = Database('sqlite:///:memory:', LRUCache, cache_size=3)
    d.create_tables()
    return d


@pytest.fixture
def database_using_lfucache():
    d = Database('sqlite:///:memory:', LFUCache, cache_size=3)
    d.create_tables()
    return d


@pytest.fixture
def database():
    d = Database('sqlite:///:memory:')
    d.create_tables()
    return d


@pytest.fixture
def query_result():
    return vso.VSOClient().query(
        vso.attrs.Time('20130801T200000', '20130801T200030'),
        vso.attrs.Instrument('PLASTIC'))


@pytest.fixture
def filled_database():
    database = Database('sqlite:///:memory:')
    database.create_tables()
    for i in xrange(1, 11):
        entry = DatabaseEntry()
        database.add(entry)
        # every fourth entry gets the tag 'foo'
        if i % 4 == 0:
            database.tag(entry, 'foo')
        # every fifth entry gets the tag 'bar'
        if i % 5 == 0:
            database.tag(entry, 'bar')
    database.commit()
    return database


def test_tags_unique(database):
    entry = DatabaseEntry()
    entry.tags = [Tag('foo')]
    database.add(entry)
    database.commit()
    entry.tags.append(Tag('foo'))
    with pytest.raises(sqlalchemy.orm.exc.FlushError):
        database.commit()


def test_setting_cache_size(database_using_lrucache):
    assert database_using_lrucache.cache_maxsize == 3
    assert database_using_lrucache.cache_size == 0
    for _ in xrange(5):
        database_using_lrucache.add(DatabaseEntry())
    assert len(database_using_lrucache) == 3
    assert database_using_lrucache.cache_size == 3
    assert database_using_lrucache.cache_maxsize == 3
    database_using_lrucache.set_cache_size(5)
    assert database_using_lrucache.cache_size == 3
    assert database_using_lrucache.cache_maxsize == 5
    for _ in xrange(5):
        database_using_lrucache.add(DatabaseEntry())
    assert len(database_using_lrucache) == 5
    assert database_using_lrucache.cache_size == 5
    assert database_using_lrucache.cache_maxsize == 5


def test_setting_cache_size_shrinking(database_using_lrucache):
    assert database_using_lrucache.cache_maxsize == 3
    assert database_using_lrucache.cache_size == 0
    for _ in xrange(5):
        database_using_lrucache.add(DatabaseEntry())
    assert len(database_using_lrucache) == 3
    assert database_using_lrucache.cache_maxsize == 3
    assert database_using_lrucache.cache_size == 3
    database_using_lrucache.set_cache_size(2)
    assert database_using_lrucache.cache_maxsize == 2
    assert database_using_lrucache.cache_size == 2
    assert len(database_using_lrucache) == 2
    assert list(database_using_lrucache) == [
        DatabaseEntry(id=4),
        DatabaseEntry(id=5)]
    for _ in xrange(5):
        database_using_lrucache.add(DatabaseEntry())
    assert len(database_using_lrucache) == 2
    assert database_using_lrucache.cache_maxsize == 2
    assert database_using_lrucache.cache_size == 2


def test_setting_cache_size_undo(database_using_lrucache):
    assert database_using_lrucache.cache_maxsize == 3
    assert database_using_lrucache.cache_size == 0
    for _ in xrange(5):
        database_using_lrucache.add(DatabaseEntry())
    assert len(database_using_lrucache) == 3
    database_using_lrucache.set_cache_size(1)
    assert database_using_lrucache.cache_size == 1
    assert len(database_using_lrucache) == 1
    database_using_lrucache.undo()
    assert len(database_using_lrucache) == 3


def test_create_tables(database_without_tables):
    assert not database_without_tables._engine.has_table('data')
    database_without_tables.create_tables()
    assert database_without_tables._engine.has_table('data')


def test_get_entry_by_id_invalid(database):
    with pytest.raises(EntryNotFoundError):
        database.get_entry_by_id(0)


def test_get_entry_by_id_zero(filled_database):
    with pytest.raises(EntryNotFoundError):
        filled_database.get_entry_by_id(0)


def test_get_entry_by_id_accessible(filled_database):
    assert filled_database.get_entry_by_id(1) == DatabaseEntry(id=1)


def test_tags_property(database):
    assert database.tags == []


def test_get_existing_tag(database):
    entry = DatabaseEntry()
    database.tag(entry, 'tag')
    database.add(entry)
    expected_tag = Tag('tag')
    expected_tag.id = 1
    assert database.get_tag('tag') == expected_tag


def test_get_nonexting_tag(database):
    with pytest.raises(NoSuchTagError):
        database.get_tag('foo')


def test_tag_missing_tags_arg(database):
    with pytest.raises(TypeError):
        database.tag(DatabaseEntry())


def test_tag_new_tag(database):
    entry = DatabaseEntry()
    database.tag(entry, 'tag')
    assert len(entry.tags) == 1
    database.add(entry)
    assert len(database.tags) == 1
    tag = entry.tags[0]
    assert tag.name == 'tag'
    assert tag in database.tags


def test_tag_existing_tag(database):
    entry1 = DatabaseEntry()
    entry2 = DatabaseEntry()
    database.tag(entry1, 'tag')
    database.tag(entry2, 'tag')
    assert entry1.tags == entry2.tags


def test_tag_duplicate(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.tag(entry, 'tag')
    database.commit()
    with pytest.raises(TagAlreadyAssignedError):
        database.tag(entry, 'tag')


def test_tag_duplicates_before_adding(database):
    entry1 = DatabaseEntry()
    entry2 = DatabaseEntry()
    database.tag(entry1, 'tag')
    database.tag(entry2, 'tag')
    database.add(entry1)
    database.add(entry2)
    with pytest.raises(sqlalchemy.orm.exc.FlushError):
        database.commit()


def remove_nonexisting_tag(database):
    with pytest.raises(NoSuchTagError):
        database.remove_tag('foo')


def test_remove_tag(filled_database):
    foo = Tag('foo')
    foo.id = 1
    fourth_entry = filled_database.get_entry_by_id(4)
    assert foo in fourth_entry.tags
    eighth_entry = filled_database.get_entry_by_id(8)
    assert foo in eighth_entry.tags
    filled_database.remove_tag(fourth_entry, 'foo')
    assert foo not in fourth_entry.tags
    assert foo in filled_database.tags
    filled_database.remove_tag(eighth_entry, 'foo')
    assert foo not in eighth_entry.tags
    assert foo not in filled_database.tags


def test_remove_tag_undo(filled_database):
    foo = Tag('foo')
    foo.id = 1
    fourth_entry = filled_database.get_entry_by_id(4)
    assert foo in fourth_entry.tags
    filled_database.remove_tag(fourth_entry, 'foo')
    assert foo not in fourth_entry.tags
    filled_database.undo()
    assert foo in fourth_entry.tags
    filled_database.redo()
    assert foo not in fourth_entry.tags
    eighth_entry = filled_database.get_entry_by_id(8)
    filled_database.remove_tag(eighth_entry, 'foo')
    assert foo not in eighth_entry.tags
    assert foo not in filled_database.tags
    filled_database.undo(2)
    assert foo not in fourth_entry.tags
    assert foo in eighth_entry.tags
    assert foo in filled_database.tags


def test_star_entry(database):
    entry = DatabaseEntry()
    assert not entry.starred
    database.star(entry)
    assert entry.starred


def test_star_already_starred_entry(database):
    entry = DatabaseEntry()
    database.star(entry)
    with pytest.raises(EntryAlreadyStarredError):
        database.star(entry)


def unstar_entry(database):
    entry = DatabaseEntry()
    assert not entry.starred
    database.star(entry)
    assert entry.starred
    database.unstar(entry)
    assert not entry.starred


def test_unstar_already_unstarred_entry(database):
    with pytest.raises(EntryAlreadyUnstarredError):
        database.unstar(DatabaseEntry())


def test_unstar_already_unstarred_entry_ignore(database):
    entry = DatabaseEntry()
    database.unstar(entry, True)
    assert not entry.starred


def test_add_entry(database):
    entry = DatabaseEntry()
    assert entry.id is None
    database.add(entry)
    database.commit()
    assert entry.id == 1


def test_add_already_existing_entry(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.commit()
    with pytest.raises(EntryAlreadyAddedError):
        database.add(entry)


def test_add_already_existing_entry_ignore(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.add(entry, True)
    database.commit()
    assert entry.id == 1


def test_add_entry_from_qr(database, query_result):
    assert len(database) == 0
    database.add_from_vso_query_result(query_result)
    assert len(database) == 4
    expected_entries = [
        DatabaseEntry(
            id=1, source="STEREO_A", provider="SSC",
            fileid="plastic/level1/ahead/2013/STA_L1_PLA_20130801_213",
            observation_time_start=datetime(2013, 8, 1, 0, 0, 0),
            observation_time_end=datetime(2013, 8, 2, 0, 0, 0),
            instrument="PLASTIC", size=166.362, waveunit="keV", wavemin=0.2,
            wavemax=100),
        DatabaseEntry(
            id=2, source="STEREO_A", provider="SSC",
            fileid="plastic/level1/ahead/2013/STA_L1_PLA_SC_20130801_213",
            observation_time_start=datetime(2013, 8, 1, 0, 0, 0),
            observation_time_end=datetime(2013, 8, 2, 0, 0, 0),
            instrument="PLASTIC", size=21.167, waveunit="keV", wavemin=0.2,
            wavemax=100),
        DatabaseEntry(
            id=3, source="STEREO_B", provider="SSC",
            fileid="plastic/level1/behind/2013/STB_L1_PLA_20130801_213",
            observation_time_start=datetime(2013, 8, 1, 0, 0, 0),
            observation_time_end=datetime(2013, 8, 2, 0, 0, 0),
            instrument="PLASTIC", size=13164.4, waveunit="keV", wavemin=0.2,
            wavemax=100),
        DatabaseEntry(
            id=4, source="STEREO_B", provider="SSC",
            fileid="plastic/level1/behind/2013/STB_L1_PLA_SC_20130801_213",
            observation_time_start=datetime(2013, 8, 1, 0, 0, 0),
            observation_time_end=datetime(2013, 8, 2, 0, 0, 0),
            instrument="PLASTIC", size=77.9795, waveunit="keV", wavemin=0.2,
            wavemax=100)]
    assert list(database) == expected_entries
    database.undo()
    assert len(database) == 0
    database.redo()
    assert len(database) == 4
    assert list(database) == expected_entries


def test_add_entries_from_qr_duplicates(database, query_result):
    assert len(database) == 0
    database.add_from_vso_query_result(query_result)
    assert len(database) == 4
    with pytest.raises(EntryAlreadyAddedError):
        database.add_from_vso_query_result(query_result)


def test_add_entries_from_qr_ignore_duplicates(database, query_result):
    assert len(database) == 0
    database.add_from_vso_query_result(query_result)
    assert len(database) == 4
    database.add_from_vso_query_result(query_result, True)
    assert len(database) == 8


def test_add_fom_path(database):
    assert len(database) == 0
    database.add_from_path(sunpy.data.sample.rootdir)
    assert len(database) == 5
    expected_entries = [
        DatabaseEntry(
            id=1, observation_time_start=datetime(2010, 10, 16, 19, 12, 18),
            observation_time_end=datetime(2010, 10, 16, 19, 12, 22),
            instrument='RHESSI',
            fits_header_entries=[
                FitsHeaderEntry('SIMPLE', True),
                FitsHeaderEntry('BITPIX', -32),
                FitsHeaderEntry('NAXIS', 2),
                FitsHeaderEntry('NAXIS1', 64),
                FitsHeaderEntry('NAXIS2', 64),
                FitsHeaderEntry('EXTEND', True),
                FitsHeaderEntry('DATE', '2011-08-29T09:50:24'),
                FitsHeaderEntry('ORIGIN', 'RHESSI'),
                FitsHeaderEntry('OBSERVER', 'schriste'),
                FitsHeaderEntry('TELESCOP', 'RHESSI'),
                FitsHeaderEntry('INSTRUME', 'RHESSI'),
                FitsHeaderEntry('OBJECT', 'Sun'),
                FitsHeaderEntry('DATE_OBS', '2010-10-16T19:12:18.000'),
                FitsHeaderEntry('DATE_END', '2010-10-16T19:12:22.000'),
                FitsHeaderEntry('TIME_UNI', 1),
                FitsHeaderEntry('ENERGY_L', 12.0),
                FitsHeaderEntry('ENERGY_H', 25.0),
                FitsHeaderEntry('TIMESYS', '1979.00'),
                FitsHeaderEntry('TIMEUNIT', 's'),
                FitsHeaderEntry('CRPIX1', -66.2278),
                FitsHeaderEntry('CRPIX2', 131.958),
                FitsHeaderEntry('CRVAL1', 0.0),
                FitsHeaderEntry('CRVAL2', 0.0),
                FitsHeaderEntry('CDELT1', 4.0),
                FitsHeaderEntry('CDELT2', 4.0),
                FitsHeaderEntry('CTYPE1', 'arcsec'),
                FitsHeaderEntry('CTYPE2', 'arcsec'),
                FitsHeaderEntry('XCEN', 394.911),
                FitsHeaderEntry('YCEN', -397.831),
                FitsHeaderEntry('CROTACN1', 0.0),
                FitsHeaderEntry('CROTACN2', 0.0),
                FitsHeaderEntry('CROTA', 0.0),
                FitsHeaderEntry('COMMENT', ''),
                FitsHeaderEntry('KEYCOMMENTS', ('  SIMPLE  Written by IDL:  Mo'
                    'n Aug 29 09:50:24 2011\n  BITPIX  Real*4 (floating point)'
                    '\n   NAXIS  \n  NAXIS1  \n  NAXIS2  \n  EXTEND  File cont'
                    'ains extensions\n    DATE  File creation date (YYYY-MM-DD'
                    'Thh:mm:ss UTC)\n  ORIGIN  High Energy Solar Spectroscopic'
                    ' Imager\nOBSERVER  Usually the name of the user who gener'
                    'ated file\nTELESCOP  Name of the Telescope or Mission\nIN'
                    'STRUME  Name of the instrument\n  OBJECT  Object being ob'
                    'served\nDATE_OBS  nominal U.T. date when integration of t'
                    'his\nDATE_END  nominal U.T. date when integration of this'
                    '\nTIME_UNI  \nENERGY_L  \nENERGY_H  \n TIMESYS  Reference'
                    ' Time\nTIMEUNIT  Does not take leap seconds into account'
                    '\n  CRPIX1  Reference pixel coordinates\n  CRPIX2  Refere'
                    'nce pixel coordinates\n  CRVAL1  Reference data coordinat'
                    'es\n  CRVAL2  Reference data coordinates\n  CDELT1  Width'
                    ' of a pixel in data units\n  CDELT2  Height of a pixel in'
                    ' data units\n  CTYPE1  data units for CDELT1\n  CTYPE2  d'
                    'ata units for CDELT2\n    XCEN  Center of image rel to su'
                    'n center ,+=W\n    YCEN  Center of image rel to sun cente'
                    'r, +=N\nCROTACN1  X Position of Center of Rotation (arcse'
                    'c)\nCROTACN2  Y Position of Center of Rotation (arcsec)\n'
                    '   CROTA  Rotation Angle (clockwise from N)')),
                FitsHeaderEntry('HISTORY', '')]),
        DatabaseEntry(
            id=2,
            observation_time_start=datetime(2002, 6, 25, 10, 0, 10, 514000),
            instrument='EIT', wavemin=195.0, wavemax=195, fits_header_entries=[
                FitsHeaderEntry('SIMPLE', True),
                FitsHeaderEntry('BITPIX', -64),
                FitsHeaderEntry('NAXIS', 2),
                FitsHeaderEntry('NAXIS1', 1024),
                FitsHeaderEntry('NAXIS2', 1024),
                FitsHeaderEntry('', ('\n\n\n\n                       / 284 = F'
                    'e XV, 304 = He II\n\n\n\n\n\n\n\n')),
                FitsHeaderEntry('DATE', '2002-06-25'),
                FitsHeaderEntry('TIME-OBS', '10:00:10'),
                FitsHeaderEntry('DATE-OBS', '2002-06-25T10:00:10.514'),
                FitsHeaderEntry('ORIGIN', 'Rocket Science'),
                FitsHeaderEntry('DATASRC', 'LZ file'),
                FitsHeaderEntry('TELESCOP', 'SOHO'),
                FitsHeaderEntry('INSTRUME', 'EIT'),
                FitsHeaderEntry('OBJECT', 'full FOV'),
                FitsHeaderEntry('BSCALE', 1.0),
                FitsHeaderEntry('BZERO', 0.0),
                FitsHeaderEntry('BUNIT', 'counts / pixel'),
                FitsHeaderEntry('WAVELNTH', 195),
                FitsHeaderEntry('FILTER', 'Al +1'),
                FitsHeaderEntry('DATE_OBS', '2002-06-25T10:00:10.514Z'),
                FitsHeaderEntry('SCI_OBJ', 'CME WATCH 195'),
                FitsHeaderEntry('OBS_PROG', '195_10S_AL_1.000'),
                FitsHeaderEntry('CMP_NO', 1),
                FitsHeaderEntry('EXPTIME', 13.298),
                FitsHeaderEntry('EXPMODE', 'backside'),
                FitsHeaderEntry('FILENAME', 'efz20020625.100010'),
                FitsHeaderEntry('CFTEMP', 0.0),
                FitsHeaderEntry('CCDTEMP', 0.0),
                FitsHeaderEntry('CTYPE1', 'Solar-X'),
                FitsHeaderEntry('CTYPE2', 'Solar-Y'),
                FitsHeaderEntry('CRPIX1', 507.01),
                FitsHeaderEntry('CRPIX2', 517.84),
                FitsHeaderEntry('CRVAL1', 0.0),
                FitsHeaderEntry('CRVAL2', 0.0),
                FitsHeaderEntry('CDELT1', 2.62),
                FitsHeaderEntry('CDELT2', 2.62),
                FitsHeaderEntry('SOLAR_R', 363.47),
                FitsHeaderEntry('SOLAR_B0', 2.17),
                FitsHeaderEntry('SC_X0', 0.0),
                FitsHeaderEntry('SC_Y0', -198.01),
                FitsHeaderEntry('SC_ROLL', 0.06),
                FitsHeaderEntry('HEC_X', 9022687.0),
                FitsHeaderEntry('HEC_Y', -150483280.0),
                FitsHeaderEntry('HEC_Z', 74309.38),
                FitsHeaderEntry('CAR_ROT', 1991.0),
                FitsHeaderEntry('COMMENT', ("CORRECTED DATE_OBS = '2002-06-25T"
                    "09:59:59.398Z'  COMMANDED EXPOSURE TIME =   10.000 s  SHU"
                    "TTER CLOSE TIME =    3.298 s  LINE_SYNC = 'no'  CAMERA_ER"
                    "R = 'yes'  IMAGE_OF_SEQ =                    0  READOUT P"
                    "ORT = 'B'  NUM_LEB_PROC =                    3  LEB_PROC "
                    "= 26 (no ROI)  LEB_PROC = 27 (no occ mask)  LEB_PROC = 12"
                    " (Rice)  BLOCKS_HORZ =   32  BLOCKS_VERT =   32  P1_X =  "
                    "         1  P2_X =        1024  P1_Y =          20  P2_Y "
                    "=        1043  N_MISSING_BLOCKS =    0")),
                FitsHeaderEntry('HISTORY', ('Version 3.3, 1999 January 5Replac'
                    'ed missing blocks with EIT_DARKSubtracted dark current, E'
                    'IT_DARKDegridded image with EIT_DEGRIDFlat fielded image '
                    'with EIT_FLATExposure normalized (per sec)Flux normalized'
                    ' to CLEAR FilterFlux normalized for time response changes'
                    )),
                FitsHeaderEntry('KEYCOMMENTS', ('  SIMPLE  Written by IDL:  Mo'
                    'n Sep  5 20:43:48 2011\n  BITPIX  IEEE double precision f'
                    'loating point\n   NAXIS  \n  NAXIS1  Number of columns\n '
                    ' NAXIS2  Number of rows\n          \n    DATE  Date of fi'
                    'le creation\nTIME-OBS  \nDATE-OBS  UTC at spacecraft\n   '
                    '       \n  ORIGIN  Rocket Science = NASA GSFC\n DATASRC  '
                    '\nTELESCOP  \nINSTRUME  \n  OBJECT  \n          \n  BSCAL'
                    'E  \n   BZERO  \n   BUNIT  \n          \nWAVELNTH  171 = '
                    'Fe IX/X, 195 = Fe XII,\n          \n  FILTER  \n         '
                    ' \nDATE_OBS  \n SCI_OBJ  \nOBS_PROG  \n  CMP_NO  Unique c'
                    'ampaign instance (1 = synoptic)\n EXPTIME  s (total comma'
                    'nded + shutter close)\n EXPMODE  \nFILENAME  \n          '
                    '\n  CFTEMP  CCD cold finger temperature (C)\n CCDTEMP  CC'
                    'D temperature (DN/100)\n          \n  CTYPE1  \n  CTYPE2 '
                    ' \n  CRPIX1  Sun center x, EIT pixels\n  CRPIX2  Sun cent'
                    'er y, EIT pixels\n  CRVAL1  \n  CRVAL2  \n  CDELT1  Pixel'
                    ' scale x (arc sec, fixed)\n  CDELT2  Pixel scale y (arc s'
                    'ec, fixed)\n SOLAR_R  Solar photospheric radius, EIT pixe'
                    'ls\nSOLAR_B0  Degrees\n          \n   SC_X0  s/c yaw (arc'
                    ' sec)\n   SC_Y0  s/c pitch (arc sec; -109.14 after 1996/0'
                    '4/16)\n SC_ROLL  s/c roll (deg., Solar north + CCW from n'
                    'ominal)\n   HEC_X  s/c heliocentric ecliptic x (km)\n   H'
                    'EC_Y  s/c heliocentric ecliptic y (km)\n   HEC_Z  s/c hel'
                    'iocentric ecliptic z (km)\n CAR_ROT  Carrington rotation '
                    'at earth\n          \n COMMENT  \n COMMENT  \n COMMENT  '
                    '\n COMMENT  \n COMMENT  \n COMMENT  \n COMMENT  \n COMMEN'
                    'T  \n COMMENT  \n COMMENT  \n COMMENT  \n COMMENT  \n COM'
                    'MENT  \n COMMENT  \n COMMENT  \n COMMENT  \n COMMENT  \n '
                    'COMMENT  \n          \n          \n HISTORY  \n          '
                    '\n HISTORY  \n HISTORY  \n HISTORY  \n HISTORY  \n HISTOR'
                    'Y  \n HISTORY  \n HISTORY  '))]),
        DatabaseEntry(
            id=3,
            observation_time_start=datetime(2011, 3, 19, 10, 54, 0, 340000),
            instrument='AIA_3', wavemin=171, wavemax=171, fits_header_entries=[
                FitsHeaderEntry('SIMPLE', True),
                FitsHeaderEntry('BITPIX', 32),
                FitsHeaderEntry('NAXIS', 2),
                FitsHeaderEntry('NAXIS1', 1024),
                FitsHeaderEntry('NAXIS2', 1024),
                FitsHeaderEntry('EXTEND', True),
                FitsHeaderEntry('COMMENT', ("FITS (Flexible Image Transport Sy"
                    "stem) format is defined in 'Astronomy  and Astrophysics',"
                    " volume 376, page 359; bibcode: 2001A&A...376..359H")),
                FitsHeaderEntry('ORIGIN', 'SDO/JSOC-SDP'),
                FitsHeaderEntry('DATE', '2011-03-19T11:08:25'),
                FitsHeaderEntry('TELESCOP', 'SDO/AIA'),
                FitsHeaderEntry('INSTRUME', 'AIA_3'),
                FitsHeaderEntry('DATE-OBS', '2011-03-19T10:54:00.34'),
                FitsHeaderEntry('T_OBS', '2011-03-19T10:54:01.34Z'),
                FitsHeaderEntry('TOBSSTEP', 90.0),
                FitsHeaderEntry('TOBSEPOC', '1977.01.01_00:00:00_TAI'),
                FitsHeaderEntry('CAMERA', 3),
                FitsHeaderEntry('IMG_TYPE', 'LIGHT'),
                FitsHeaderEntry('EXPTIME', 1.999601),
                FitsHeaderEntry('EXPSDEV', 0.00016),
                FitsHeaderEntry('INT_TIME', 2.273438),
                FitsHeaderEntry('WAVELNTH', 171),
                FitsHeaderEntry('WAVEUNIT', 'angstrom'),
                FitsHeaderEntry('WAVE_STR', '171_THIN'),
                FitsHeaderEntry('FSN', 22642033),
                FitsHeaderEntry('FID', 0),
                FitsHeaderEntry('LVL_NUM', 1.5),
                FitsHeaderEntry('QUALLEV0', 0),
                FitsHeaderEntry('QUALITY', 1073741824),
                FitsHeaderEntry('TOTVALS', 16777216),
                FitsHeaderEntry('DATAVALS', 16777216),
                FitsHeaderEntry('MISSVALS', 0),
                FitsHeaderEntry('PERCENTD', 100.0),
                FitsHeaderEntry('DATAMIN', -7),
                FitsHeaderEntry('DATAMAX', 16228),
                FitsHeaderEntry('DATAMEDN', 192),
                FitsHeaderEntry('DATAMEAN', 236.57),
                FitsHeaderEntry('DATARMS', 294.02),
                FitsHeaderEntry('DATASKEW', 4.63),
                FitsHeaderEntry('DATAKURT', 49.56),
                FitsHeaderEntry('OSCNMEAN', 'nan'),
                FitsHeaderEntry('OSCNRMS', 'nan'),
                FitsHeaderEntry('FLAT_REC', 'aia.flatfield[:#7]'),
                FitsHeaderEntry('CTYPE1', 'HPLN-TAN'),
                FitsHeaderEntry('CUNIT1', 'arcsec'),
                FitsHeaderEntry('CRVAL1', 0.0),
                FitsHeaderEntry('CDELT1', 2.4),
                FitsHeaderEntry('CRPIX1', 512.5),
                FitsHeaderEntry('CTYPE2', 'HPLT-TAN'),
                FitsHeaderEntry('CUNIT2', 'arcsec'),
                FitsHeaderEntry('CRVAL2', 0.0),
                FitsHeaderEntry('CDELT2', 2.4),
                FitsHeaderEntry('CRPIX2', 512.5),
                FitsHeaderEntry('CROTA2', 0.0),
                FitsHeaderEntry('R_SUN', 1608.94397),
                FitsHeaderEntry('MPO_REC', 'sdo.master_pointing[:#116]'),
                FitsHeaderEntry('INST_ROT', 0.102488),
                FitsHeaderEntry('IMSCL_MP', 0.599076),
                FitsHeaderEntry('X0_MP', 2052.399902),
                FitsHeaderEntry('Y0_MP', 2046.589966),
                FitsHeaderEntry('RSUN_LF', 'nan'),
                FitsHeaderEntry('X0_LF', 'nan'),
                FitsHeaderEntry('Y0_LF', 'nan'),
                FitsHeaderEntry('ASD_REC', 'sdo.lev0_asd_0004[:#8948067]'),
                FitsHeaderEntry('SAT_Y0', -0.365593),
                FitsHeaderEntry('SAT_Z0', 14.820671),
                FitsHeaderEntry('SAT_ROT', -8.8e-05),
                FitsHeaderEntry('ACS_MODE', 'SCIENCE'),
                FitsHeaderEntry('ACS_ECLP', 'NO'),
                FitsHeaderEntry('ACS_SUNP', 'YES'),
                FitsHeaderEntry('ACS_SAFE', 'NO'),
                FitsHeaderEntry('ACS_CGT', 'GT3'),
                FitsHeaderEntry('ORB_REC',
                    'sdo.fds_orbit_vectors[2011.03.19_10:54:00_UTC]'),
                FitsHeaderEntry('DSUN_REF', 149597870691.0),
                FitsHeaderEntry('DSUN_OBS', 148940609626.98),
                FitsHeaderEntry('RSUN_REF', 696000000.0),
                FitsHeaderEntry('RSUN_OBS', 963.879683),
                FitsHeaderEntry('GCIEC_X', 'nan'),
                FitsHeaderEntry('GCIEC_Y', 'nan'),
                FitsHeaderEntry('GCIEC_Z', 'nan'),
                FitsHeaderEntry('HCIEC_X', 'nan'),
                FitsHeaderEntry('HCIEC_Y', 'nan'),
                FitsHeaderEntry('HCIEC_Z', 'nan'),
                FitsHeaderEntry('OBS_VR', -2132.568376),
                FitsHeaderEntry('OBS_VW', 30775.731671),
                FitsHeaderEntry('OBS_VN', 2177.6711),
                FitsHeaderEntry('CRLN_OBS', 315.285065),
                FitsHeaderEntry('CRLT_OBS', -7.064078),
                FitsHeaderEntry('CAR_ROT', 2108),
                FitsHeaderEntry('ROI_NWIN', -2147483648),
                FitsHeaderEntry('ROI_SUM', -2147483648),
                FitsHeaderEntry('ROI_NAX1', -2147483648),
                FitsHeaderEntry('ROI_NAY1', -2147483648),
                FitsHeaderEntry('ROI_LLX1', -2147483648),
                FitsHeaderEntry('ROI_LLY1', -2147483648),
                FitsHeaderEntry('ROI_NAX2', -2147483648),
                FitsHeaderEntry('ROI_NAY2', -2147483648),
                FitsHeaderEntry('ROI_LLX2', -2147483648),
                FitsHeaderEntry('ROI_LLY2', -2147483648),
                FitsHeaderEntry('ISPSNAME', 'aia.lev0_isp_0011'),
                FitsHeaderEntry('ISPPKTIM', '2011-03-19T10:53:57.50Z'),
                FitsHeaderEntry('ISPPKTVN', '001.197'),
                FitsHeaderEntry('AIVNMST', 453),
                FitsHeaderEntry('AIMGOTS', 1679223275),
                FitsHeaderEntry('ASQHDR', 2170125681L),
                FitsHeaderEntry('ASQTNUM', 2),
                FitsHeaderEntry('ASQFSN', 22642033),
                FitsHeaderEntry('AIAHFSN', 22642025),
                FitsHeaderEntry('AECDELAY', 1540),
                FitsHeaderEntry('AIAECTI', 0),
                FitsHeaderEntry('AIASEN', 0),
                FitsHeaderEntry('AIFDBID', 241),
                FitsHeaderEntry('AIMGOTSS', 5339),
                FitsHeaderEntry('AIFCPS', 10),
                FitsHeaderEntry('AIFTSWTH', 0),
                FitsHeaderEntry('AIFRMLID', 3024),
                FitsHeaderEntry('AIFTSID', 40960),
                FitsHeaderEntry('AIHISMXB', 7),
                FitsHeaderEntry('AIHIS192', 8381297),
                FitsHeaderEntry('AIHIS348', 8388262),
                FitsHeaderEntry('AIHIS604', 8388597),
                FitsHeaderEntry('AIHIS860', 8388603),
                FitsHeaderEntry('AIFWEN', 204),
                FitsHeaderEntry('AIMGSHCE', 2000),
                FitsHeaderEntry('AECTYPE', 0),
                FitsHeaderEntry('AECMODE', 'ON'),
                FitsHeaderEntry('AISTATE', 'CLOSED'),
                FitsHeaderEntry('AIAECENF', 1),
                FitsHeaderEntry('AIFILTYP', 0),
                FitsHeaderEntry('AIMSHOBC', 41.099998),
                FitsHeaderEntry('AIMSHOBE', 26.076),
                FitsHeaderEntry('AIMSHOTC', 55.312),
                FitsHeaderEntry('AIMSHOTE', 69.316002),
                FitsHeaderEntry('AIMSHCBC', 2040.791992),
                FitsHeaderEntry('AIMSHCBE', 2025.864014),
                FitsHeaderEntry('AIMSHCTC', 2054.875977),
                FitsHeaderEntry('AIMSHCTE', 2068.676025),
                FitsHeaderEntry('AICFGDL1', 0),
                FitsHeaderEntry('AICFGDL2', 137),
                FitsHeaderEntry('AICFGDL3', 201),
                FitsHeaderEntry('AICFGDL4', 236),
                FitsHeaderEntry('AIFOENFL', 1),
                FitsHeaderEntry('AIMGFSN', 5),
                FitsHeaderEntry('AIMGTYP', 0),
                FitsHeaderEntry('AIAWVLEN', 7),
                FitsHeaderEntry('AIAGP1', 0),
                FitsHeaderEntry('AIAGP2', 0),
                FitsHeaderEntry('AIAGP3', 0),
                FitsHeaderEntry('AIAGP4', 0),
                FitsHeaderEntry('AIAGP5', 0),
                FitsHeaderEntry('AIAGP6', 0),
                FitsHeaderEntry('AIAGP7', 0),
                FitsHeaderEntry('AIAGP8', 393),
                FitsHeaderEntry('AIAGP9', 457),
                FitsHeaderEntry('AIAGP10', 748),
                FitsHeaderEntry('AGT1SVY', 18),
                FitsHeaderEntry('AGT1SVZ', 10),
                FitsHeaderEntry('AGT2SVY', 10),
                FitsHeaderEntry('AGT2SVZ', 8),
                FitsHeaderEntry('AGT3SVY', 2),
                FitsHeaderEntry('AGT3SVZ', 0),
                FitsHeaderEntry('AGT4SVY', 14),
                FitsHeaderEntry('AGT4SVZ', 0),
                FitsHeaderEntry('AIMGSHEN', 4),
                FitsHeaderEntry('RECNUM', 76202),
                FitsHeaderEntry('BLANK', -2147483648),
                FitsHeaderEntry('BZERO', 0.0),
                FitsHeaderEntry('BSCALE', 0.0625),
                FitsHeaderEntry('CHECKSUM', 'J7qAL7o6J7oAJ7o5'),
                FitsHeaderEntry('DATASUM', '3958014355'),
                FitsHeaderEntry('KEYCOMMENTS', ('  SIMPLE  file does conform '
                    'to FITS standard\n  BITPIX  number of bits per data pixel'
                    '\n   NAXIS  number of data axes\n  NAXIS1  length of data'
                    ' axis 1\n  NAXIS2  length of data axis 2\n  EXTEND  FITS '
                    'dataset may contain extensions\n COMMENT  \n COMMENT  \n '
                    ' ORIGIN  \n    DATE  \nTELESCOP  \nINSTRUME  \nDATE-OBS  '
                    '\n   T_OBS  \nTOBSSTEP  \nTOBSEPOC  \n  CAMERA  \nIMG_TYP'
                    'E  \n EXPTIME  \n EXPSDEV  \nINT_TIME  \nWAVELNTH  \nWAVE'
                    'UNIT  \nWAVE_STR  \n     FSN  \n     FID  \n LVL_NUM  \nQ'
                    'UALLEV0  \n QUALITY  \n TOTVALS  \nDATAVALS  \nMISSVALS  '
                    '\nPERCENTD  \n DATAMIN  \n DATAMAX  \nDATAMEDN  \nDATAMEA'
                    'N  \n DATARMS  \nDATASKEW  \nDATAKURT  \nOSCNMEAN  \n OSC'
                    'NRMS  \nFLAT_REC  \n  CTYPE1  \n  CUNIT1  \n  CRVAL1  \n '
                    ' CDELT1  \n  CRPIX1  \n  CTYPE2  \n  CUNIT2  \n  CRVAL2  '
                    '\n  CDELT2  \n  CRPIX2  \n  CROTA2  \n   R_SUN  \n MPO_RE'
                    'C  \nINST_ROT  \nIMSCL_MP  \n   X0_MP  \n   Y0_MP  \n RSU'
                    'N_LF  \n   X0_LF  \n   Y0_LF  \n ASD_REC  \n  SAT_Y0  \n '
                    ' SAT_Z0  \n SAT_ROT  \nACS_MODE  \nACS_ECLP  \nACS_SUNP  '
                    '\nACS_SAFE  \n ACS_CGT  \n ORB_REC  \nDSUN_REF  \nDSUN_OB'
                    'S  \nRSUN_REF  \nRSUN_OBS  \n GCIEC_X  \n GCIEC_Y  \n GCI'
                    'EC_Z  \n HCIEC_X  \n HCIEC_Y  \n HCIEC_Z  \n  OBS_VR  \n '
                    ' OBS_VW  \n  OBS_VN  \nCRLN_OBS  \nCRLT_OBS  \n CAR_ROT  '
                    '\nROI_NWIN  \n ROI_SUM  \nROI_NAX1  \nROI_NAY1  \nROI_LLX'
                    '1  \nROI_LLY1  \nROI_NAX2  \nROI_NAY2  \nROI_LLX2  \nROI_'
                    'LLY2  \nISPSNAME  \nISPPKTIM  \nISPPKTVN  \n AIVNMST  \n '
                    'AIMGOTS  \n  ASQHDR  \n ASQTNUM  \n  ASQFSN  \n AIAHFSN  '
                    '\nAECDELAY  \n AIAECTI  \n  AIASEN  \n AIFDBID  \nAIMGOTS'
                    'S  \n  AIFCPS  \nAIFTSWTH  \nAIFRMLID  \n AIFTSID  \nAIHI'
                    'SMXB  \nAIHIS192  \nAIHIS348  \nAIHIS604  \nAIHIS860  \n '
                    ' AIFWEN  \nAIMGSHCE  \n AECTYPE  \n AECMODE  \n AISTATE  '
                    '\nAIAECENF  \nAIFILTYP  \nAIMSHOBC  \nAIMSHOBE  \nAIMSHOT'
                    'C  \nAIMSHOTE  \nAIMSHCBC  \nAIMSHCBE  \nAIMSHCTC  \nAIMS'
                    'HCTE  \nAICFGDL1  \nAICFGDL2  \nAICFGDL3  \nAICFGDL4  \nA'
                    'IFOENFL  \n AIMGFSN  \n AIMGTYP  \nAIAWVLEN  \n  AIAGP1  '
                    '\n  AIAGP2  \n  AIAGP3  \n  AIAGP4  \n  AIAGP5  \n  AIAGP'
                    '6  \n  AIAGP7  \n  AIAGP8  \n  AIAGP9  \n AIAGP10  \n AGT'
                    '1SVY  \n AGT1SVZ  \n AGT2SVY  \n AGT2SVZ  \n AGT3SVY  \n '
                    'AGT3SVZ  \n AGT4SVY  \n AGT4SVZ  \nAIMGSHEN  \n  RECNUM  '
                    '\n   BLANK  \n   BZERO  \n  BSCALE  \nCHECKSUM  HDU check'
                    'sum updated 2011-03-19T11:08:18\n DATASUM  data unit chec'
                    'ksum updated 2011-03-19T11:08:18')),
                FitsHeaderEntry('HISTORY', '')]),
        DatabaseEntry(
            id=4, observation_time_start=datetime(2011, 9, 22, 0, 0, 0),
            observation_time_end=datetime(2011, 9, 22, 0, 0, 0),
            instrument='BIR',
            fits_header_entries=[
                FitsHeaderEntry('SIMPLE', True),
                FitsHeaderEntry('BITPIX', 8),
                FitsHeaderEntry('NAXIS', 2),
                FitsHeaderEntry('NAXIS1', 3600),
                FitsHeaderEntry('NAXIS2', 200),
                FitsHeaderEntry('EXTEND', True),
                FitsHeaderEntry('COMMENT', ("= 'Warning: the value of CDELT1 m"
                    "ay be rounded!'= 'Warning: the frequency axis may not be "
                    "regular!'= 'Warning: the value of CDELT2 may be rounded!'"
                    "= 'FITS Definition document #100 and other FITS informati"
                    "on'")),
                FitsHeaderEntry('DATE', '2011-09-22'),
                FitsHeaderEntry('CONTENT',
                    '2011/09/22  Radio flux density, e-CALLISTO (BIR)'),
                FitsHeaderEntry('ORIGIN', 'Birr_Castle_Ireland'),
                FitsHeaderEntry('TELESCOP', 'Radio Spectrometer'),
                FitsHeaderEntry('INSTRUME', 'BIR'),
                FitsHeaderEntry('OBJECT', 'Sun'),
                FitsHeaderEntry('DATE-OBS', '2011/09/22'),
                FitsHeaderEntry('TIME-OBS', '10:30:00.051'),
                FitsHeaderEntry('DATE-END', '2011/09/22'),
                FitsHeaderEntry('TIME-END', '10:45:00'),
                FitsHeaderEntry('BZERO', 0.0),
                FitsHeaderEntry('BSCALE', 1.0),
                FitsHeaderEntry('BUNIT', 'digits'),
                FitsHeaderEntry('DATAMIN', 119),
                FitsHeaderEntry('DATAMAX', 206),
                FitsHeaderEntry('CRVAL1', 37800.0),
                FitsHeaderEntry('CRPIX1', 0),
                FitsHeaderEntry('CTYPE1', 'Time [UT]'),
                FitsHeaderEntry('CDELT1', 0.25),
                FitsHeaderEntry('CRVAL2', 200.0),
                FitsHeaderEntry('CRPIX2', 0),
                FitsHeaderEntry('CTYPE2', 'Frequency [MHz]'),
                FitsHeaderEntry('CDELT2', -1.0),
                FitsHeaderEntry('OBS_LAT', 53.0941390991211),
                FitsHeaderEntry('OBS_LAC', 'N'),
                FitsHeaderEntry('OBS_LON', 7.92012977600098),
                FitsHeaderEntry('OBS_LOC', 'E'),
                FitsHeaderEntry('OBS_ALT', 416.5),
                FitsHeaderEntry('FRQFILE', 'FRQ_TEST.CFG'),
                FitsHeaderEntry('PWM_VAL', 80),
                FitsHeaderEntry('HISTORY', "= '        '"),
                FitsHeaderEntry('KEYCOMMENTS', ('  SIMPLE  file does conform '
                    'to FITS standard\n  BITPIX  number of bits per data pixel'
                    '\n   NAXIS  number of data axes\n  NAXIS1  length of data'
                    ' axis 1\n  NAXIS2  length of data axis 2\n  EXTEND  FITS '
                    'dataset may contain extensions\n COMMENT  \n COMMENT  \n '
                    'COMMENT  \n COMMENT  \n    DATE  Time of observation\n CO'
                    'NTENT  Title of image\n  ORIGIN  Organization name\nTELES'
                    'COP  Type of instrument\nINSTRUME  Name of the spectromet'
                    'er\n  OBJECT  object description\nDATE-OBS  Date observat'
                    'ion starts\nTIME-OBS  Time observation starts\nDATE-END  '
                    'date observation ends\nTIME-END  time observation ends\n '
                    '  BZERO  scaling offset\n  BSCALE  scaling factor\n   BUN'
                    'IT  z-axis title\n DATAMIN  minimum element in image\n DA'
                    'TAMAX  maximum element in image\n  CRVAL1  value on axis '
                    '1 at reference pixel [sec of day]\n  CRPIX1  reference pi'
                    'xel of axis 1\n  CTYPE1  title of axis 1\n  CDELT1  step '
                    'between first and second element in x-axis\n  CRVAL2  val'
                    'ue on axis 2 at the reference pixel\n  CRPIX2  reference '
                    'pixel of axis 2\n  CTYPE2  title of axis 2\n  CDELT2  ste'
                    'p between first and second element in axis\n OBS_LAT  obs'
                    'ervatory latitude in degree\n OBS_LAC  observatory latitu'
                    'de code {N,S}\n OBS_LON  observatory longitude in degree'
                    '\n OBS_LOC  observatory longitude code {E,W}\n OBS_ALT  '
                    'observatory altitude in meter asl\n FRQFILE  name of freq'
                    'uency file\n PWM_VAL  PWM value to control tuner gain\n '
                    'HISTORY  '))]),
        DatabaseEntry(
            id=5, observation_time_start=datetime(2002, 2, 20, 11, 6, 0),
            observation_time_end=datetime(2002, 2, 20, 11, 6, 43, 330000),
            instrument='RHESSI',
            fits_header_entries=[
                FitsHeaderEntry('SIMPLE', True),
                FitsHeaderEntry('BITPIX', 8),
                FitsHeaderEntry('NAXIS', 0),
                FitsHeaderEntry('EXTEND', True),
                FitsHeaderEntry('DATE', '2011-09-13T15:37:38'),
                FitsHeaderEntry('ORIGIN', 'RHESSI'),
                FitsHeaderEntry('OBSERVER', 'Unknown'),
                FitsHeaderEntry('TELESCOP', 'RHESSI'),
                FitsHeaderEntry('INSTRUME', 'RHESSI'),
                FitsHeaderEntry('OBJECT', 'Sun'),
                FitsHeaderEntry('DATE_OBS', '2002-02-20T11:06:00.000'),
                FitsHeaderEntry('DATE_END', '2002-02-20T11:06:43.330'),
                FitsHeaderEntry('TIME_UNI', 1),
                FitsHeaderEntry('ENERGY_L', 25.0),
                FitsHeaderEntry('ENERGY_H', 40.0),
                FitsHeaderEntry('TIMESYS', '1979-01-01T00:00:00'),
                FitsHeaderEntry('TIMEUNIT', 'd'),
                FitsHeaderEntry('COMMENT', ''),
                FitsHeaderEntry('KEYCOMMENTS', ('  SIMPLE  Written by IDL:  Tu'
                    'e Sep 13 15:37:38 2011\n  BITPIX  \n   NAXIS  \n  EXTEND '
                    ' File contains extensions\n    DATE  File creation date ('
                    'YYYY-MM-DDThh:mm:ss UTC)\n  ORIGIN  High Energy Solar Spe'
                    'ctroscopic Imager\nOBSERVER  Usually the name of the user'
                    ' who generated file\nTELESCOP  Name of the Telescope or M'
                    'ission\nINSTRUME  Name of the instrument\n  OBJECT  Objec'
                    't being observed\nDATE_OBS  nominal U.T. date when integr'
                    'ation of this\nDATE_END  nominal U.T. date when integrati'
                    'on of this\nTIME_UNI  \nENERGY_L  \nENERGY_H  \n TIMESYS '
                    ' Reference time in YYYY MM DD hh:mm:ss\nTIMEUNIT  Unit fo'
                    'r TIMEZERO, TSTARTI and TSTOPI')),
                FitsHeaderEntry('HISTORY', '')])]
    assert list(database) == expected_entries
    database.undo()
    assert len(database) == 0
    database.redo()
    assert len(database) == 5
    assert list(database) == expected_entries


def test_add_fom_path_duplicates(database):
    database.add_from_path(sunpy.data.sample.rootdir)
    assert len(database) == 5
    with pytest.raises(EntryAlreadyAddedError):
        database.add_from_path(sunpy.data.sample.rootdir)


def test_add_fom_path_ignore_duplicates(database):
    database.add_from_path(sunpy.data.sample.rootdir)
    assert len(database) == 5
    database.add_from_path(
        sunpy.data.sample.rootdir, ignore_already_added=True)
    assert len(database) == 10


def test_edit_entry(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.commit()
    assert entry.id == 1
    database.edit(entry, id=42)
    assert entry.id == 42


def test_remove_existing_entry(database):
    entry = DatabaseEntry()
    database.add(entry)
    assert database.session.query(DatabaseEntry).count() == 1
    assert entry.id == 1
    database.remove(entry)
    assert database.session.query(DatabaseEntry).count() == 0


def test_remove_nonexisting_entry(database):
    with pytest.raises(NoSuchEntryError):
        database.remove(DatabaseEntry())


def test_clear_empty_database(database):
    database.clear()


def test_clear_database(filled_database):
    assert len(filled_database) == 10
    filled_database.clear()
    assert not filled_database
    filled_database.undo()
    assert len(filled_database) == 10
    filled_database.commit()
    filled_database.redo()
    assert not filled_database


def test_getitem_notfound(database):
    with pytest.raises(IndexError):
        database[23]


def test_getitem_one(filled_database):
    first_entry = DatabaseEntry(id=1)
    assert filled_database[0] == first_entry


def test_getitem_getall(filled_database):
    entries = filled_database[:]
    assert entries == list(filled_database)


def test_getitem_custom(filled_database):
    entries = filled_database[1:5:2]
    foo = Tag('foo')
    foo.id = 1
    assert entries == [
        DatabaseEntry(id=2), DatabaseEntry(id=4, tags=[foo])]


def test_getitem_exceeding_range(filled_database):
    entries = filled_database[7:1000]
    foo = Tag('foo')
    foo.id = 1
    bar = Tag('bar')
    bar.id = 2
    assert entries == [
        DatabaseEntry(id=8, tags=[foo]),
        DatabaseEntry(id=9),
        DatabaseEntry(id=10, tags=[bar])]


def test_contains_exists(database):
    entry = DatabaseEntry()
    database.add(entry)
    database.commit()
    assert entry in database


def test_contains_precommit(database):
    entry = DatabaseEntry()
    database.add(entry)
    assert entry not in database


def test_contains_notexists(database):
    assert DatabaseEntry() not in database


def test_iter(database):
    entry1 = DatabaseEntry()
    entry2 = DatabaseEntry()
    database.add(entry1)
    database.add(entry2)
    expected_entries = [entry1, entry2]
    entries = list(database)
    assert entries == expected_entries


def test_len(database):
    assert len(database) == 0
    database.session.add(DatabaseEntry())
    assert len(database) == 1


def test_lru_cache(database_using_lrucache):
    assert not database_using_lrucache._cache
    entry1, entry2, entry3 = DatabaseEntry(), DatabaseEntry(), DatabaseEntry()
    database_using_lrucache.add(entry1)
    database_using_lrucache.add(entry2)
    database_using_lrucache.add(entry3)
    assert len(database_using_lrucache) == 3
    entries = list(database_using_lrucache)
    assert entries[0] == entry1
    assert entries[1] == entry2
    assert entries[2] == entry3
    #assert database_using_lrucache._cache.items() == [
    #    (1, entry1), (2, entry2), (3, entry3)]
    database_using_lrucache.get_entry_by_id(1)
    database_using_lrucache.get_entry_by_id(3)
    entry4 = DatabaseEntry()
    database_using_lrucache.add(entry4)
    assert len(database_using_lrucache) == 3
    entries = list(database_using_lrucache)
    assert entries[0] == entry1
    assert entries[1] == entry3
    assert entries[2] == entry4
    #assert database_using_lrucache._cache.items() == [
    #    (1, entry1), (3, entry3), (4, entry4)]


def test_lfu_cache(database_using_lfucache):
    assert not database_using_lfucache._cache
    entry1, entry2, entry3 = DatabaseEntry(), DatabaseEntry(), DatabaseEntry()
    database_using_lfucache.add(entry1)
    database_using_lfucache.add(entry2)
    database_using_lfucache.add(entry3)
    assert len(database_using_lfucache) == 3
    entries = list(database_using_lfucache)
    assert entries[0] == entry1
    assert entries[1] == entry2
    assert entries[2] == entry3
    #assert database_using_lrucache._cache.items() == [
    #    (1, entry1), (2, entry2), (3, entry3)]
    # access the entries #1 and #2 to increment their counters
    database_using_lfucache.get_entry_by_id(1)
    database_using_lfucache.get_entry_by_id(2)
    entry4 = DatabaseEntry()
    database_using_lfucache.add(entry4)
    assert len(database_using_lfucache) == 3
    entries = list(database_using_lfucache)
    assert entries[0] == entry1
    assert entries[1] == entry2
    assert entries[2] == entry4
    #assert database_using_lrucache._cache.items() == [
    #    (1, entry1), (2, entry2), (4, entry4)]


def test_query_missing_arg(database):
    with pytest.raises(TypeError):
        database.query()


def test_query_unexpected_kwarg(database):
    with pytest.raises(TypeError):
        database.query(attrs.Starred(), foo=42)


def test_query(filled_database):
    foo = Tag('foo')
    foo.id = 1
    bar = Tag('bar')
    bar.id = 2
    entries = filled_database.query(attrs.Tag('foo') | attrs.Tag('bar'), sortby='id')
    assert len(entries) == 4
    assert entries == [
        DatabaseEntry(id=4, tags=[foo]),
        DatabaseEntry(id=5, tags=[bar]),
        DatabaseEntry(id=8, tags=[foo]),
        DatabaseEntry(id=10, tags=[bar])]
