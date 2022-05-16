import os
import pathlib
import tempfile

import numpy as np
import pytest

from astropy.io import fits
from astropy.wcs import WCS

import sunpy
import sunpy.map
from sunpy.data.test import get_dummy_map_from_header, get_test_filepath, rootdir, test_data_filenames
from sunpy.tests.helpers import skip_glymur
from sunpy.util.exceptions import NoMapsInFileError, SunpyMetadataWarning, SunpyUserWarning

a_list_of_many = [f for f in test_data_filenames() if 'efz' in f.name]

AIA_171_IMAGE = get_test_filepath('aia_171_level1.fits')
RHESSI_IMAGE = get_test_filepath('hsi_image_20101016_191218.fits')
AIA_193_JP2 = get_test_filepath("2013_06_24__17_31_30_84__SDO_AIA_AIA_193.jp2")
AIA_MAP = sunpy.map.Map(AIA_171_IMAGE)
VALID_MAP_INPUTS = [
    (AIA_171_IMAGE, ),
    (pathlib.Path(AIA_171_IMAGE), ),
    (rootdir / "EIT", ),
    (os.fspath(rootdir / "EIT"), ),
    (rootdir / "EIT" / "*.fits", ),
    (AIA_MAP, ),
    (AIA_MAP.data, AIA_MAP.meta),
    ((AIA_MAP.data, AIA_MAP.meta), ),
]


@pytest.mark.parametrize('args1', VALID_MAP_INPUTS)
@pytest.mark.parametrize('args2', VALID_MAP_INPUTS)
def test_two_map_inputs(args1, args2):
    out = sunpy.map.Map(*args1, *args2)
    if isinstance(out, list):
        for m in out:
            assert isinstance(m, sunpy.map.GenericMap)
    else:
        assert isinstance(out, sunpy.map.GenericMap)


def test_mapsequence(eit_fits_directory):
    # Test making a MapSequence
    sequence = sunpy.map.Map(list(eit_fits_directory.glob('*.fits')), sequence=True)
    assert isinstance(sequence, sunpy.map.MapSequence)


def test_mapsequence_sortby(eit_fits_directory):
    # Test making a MapSequence with sortby kwarg
    sequence = sunpy.map.Map(list(eit_fits_directory.glob('*.fits')), sequence=True, sortby=None)
    assert isinstance(sequence, sunpy.map.MapSequence)


def test_composite():
    # Test making a CompositeMap
    comp = sunpy.map.Map(AIA_171_IMAGE, RHESSI_IMAGE, composite=True)
    assert isinstance(comp, sunpy.map.CompositeMap)

# Want to check that patterns work, so ignore this warning that comes from
# the AIA test data


@pytest.mark.filterwarnings("ignore:Invalid 'BLANK' keyword in header")
def test_patterns(eit_fits_directory):
    # Test different Map pattern matching

    # File name
    aiamap = sunpy.map.Map(AIA_171_IMAGE)
    assert isinstance(aiamap, sunpy.map.GenericMap)

    # Directory
    maps = sunpy.map.Map(os.fspath(eit_fits_directory))
    assert isinstance(maps, list)
    assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])
    # Test that returned maps are sorted
    files_sorted = sorted(list(eit_fits_directory.glob('*')))
    maps_sorted = [sunpy.map.Map(os.fspath(f)) for f in files_sorted]
    assert all([m.date == m_s.date for m, m_s in zip(maps, maps_sorted)])

    # Pathlib
    path = pathlib.Path(AIA_171_IMAGE)
    aiamap = sunpy.map.Map(path)
    assert isinstance(aiamap, sunpy.map.GenericMap)
    maps = sunpy.map.Map(eit_fits_directory)
    assert isinstance(maps, list)
    assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])

    # Glob
    pattern = os.path.join(eit_fits_directory, "*")
    maps = sunpy.map.Map(pattern)
    assert isinstance(maps, list)
    assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])
    # Test that returned maps are sorted
    files_sorted = sorted(list(pathlib.Path(pattern).parent.glob('*')))
    maps_sorted = [sunpy.map.Map(os.fspath(f)) for f in files_sorted]
    assert all([m.date == m_s.date for m, m_s in zip(maps, maps_sorted)])
    # Single character wildcard (?)
    pattern = os.path.join(eit_fits_directory, "efz20040301.0?0010_s.fits")
    maps = sunpy.map.Map(pattern)
    assert isinstance(maps, list)
    assert len(maps) == 7
    assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])
    # Character ranges
    pattern = os.path.join(eit_fits_directory, "efz20040301.0[2-6]0010_s.fits")
    maps = sunpy.map.Map(pattern)
    assert isinstance(maps, list)
    assert len(maps) == 4
    assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])

    # Already a Map
    amap = sunpy.map.Map(maps[0])
    assert isinstance(amap, sunpy.map.GenericMap)

    # A list of filenames
    maps = sunpy.map.Map(list(eit_fits_directory.glob('*.fits')))
    assert isinstance(maps, list)
    assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])

    # Data-header pair in a tuple
    pair_map = sunpy.map.Map((amap.data, amap.meta))
    assert isinstance(pair_map, sunpy.map.GenericMap)

    # Data-header pair not in a tuple
    pair_map = sunpy.map.Map(amap.data, amap.meta)
    assert isinstance(pair_map, sunpy.map.GenericMap)

    # Data-wcs object pair in tuple
    pair_map = sunpy.map.Map((amap.data, WCS(AIA_171_IMAGE)))
    assert isinstance(pair_map, sunpy.map.GenericMap)

    # Data-wcs object pair not in a tuple
    pair_map = sunpy.map.Map(amap.data, WCS(AIA_171_IMAGE))
    assert isinstance(pair_map, sunpy.map.GenericMap)

    # Data-header from FITS
    with fits.open(AIA_171_IMAGE) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    pair_map = sunpy.map.Map((data, header))
    assert isinstance(pair_map, sunpy.map.GenericMap)
    pair_map, pair_map = sunpy.map.Map(((data, header), (data, header)))
    assert isinstance(pair_map, sunpy.map.GenericMap)
    pair_map = sunpy.map.Map(data, header)
    assert isinstance(pair_map, sunpy.map.GenericMap)

    # Custom Map
    data = np.arange(0, 100).reshape(10, 10)
    header = {'cdelt1': 10, 'cdelt2': 10,
              'telescop': 'sunpy',
              'cunit1': 'arcsec', 'cunit2': 'arcsec'}
    with pytest.warns(SunpyMetadataWarning, match='Missing CTYPE1 from metadata, assuming CTYPE1 is HPLN-TAN'):
        pair_map = sunpy.map.Map(data, header)
    assert isinstance(pair_map, sunpy.map.GenericMap)

    # Common keys not strings
    data = np.arange(0, 100).reshape(10, 10)
    header = {'cdelt1': 10, 'cdelt2': 10,
              'telescop': 100,
              'detector': 1,
              'instrume': 50,
              'cunit1': 'arcsec', 'cunit2': 'arcsec'}
    with pytest.warns(SunpyMetadataWarning, match='Missing CTYPE1 from metadata, assuming CTYPE1 is HPLN-TAN'):
        pair_map = sunpy.map.Map(data, header)
    assert isinstance(pair_map, sunpy.map.GenericMap)


def test_errors(tmpdir):
    # If directory doesn't exist, make sure it's listed in the error msg
    nonexist_dir = 'nonexist'
    directory = pathlib.Path(tmpdir, nonexist_dir)
    with pytest.raises(ValueError, match=nonexist_dir):
        sunpy.map.Map(os.fspath(directory))

    with pytest.raises(ValueError, match='Invalid input: 78'):
        # Check a random unsupported type (int) fails
        sunpy.map.Map(78)

    # If one file failed to load, make sure it's raised as an expection.
    p = tmpdir.mkdir("sub").join("hello.fits")
    p.write("content")
    files = [AIA_171_IMAGE, p.strpath]
    with pytest.raises(OSError, match=(fr"Failed to read *")):
        sunpy.map.Map(files)


@pytest.mark.filterwarnings("ignore:One of the data, header pairs failed to validate")
@pytest.mark.parametrize('silence,error,match',
                         [(True, RuntimeError, 'No maps loaded'),
                             (False, sunpy.map.mapbase.MapMetaValidationError,
                              'Image coordinate units for axis 1 not present in metadata.')])
def test_silence_errors(silence, error, match):
    # Check that the correct errors are raised depending on silence_errors value
    data = np.arange(0, 100).reshape(10, 10)
    header = {}
    with pytest.raises(error, match=match):
        _ = sunpy.map.Map(data, header, silence_errors=silence)


def test_dask_array():
    dask_array = pytest.importorskip('dask.array')
    amap = sunpy.map.Map(AIA_171_IMAGE)
    da = dask_array.from_array(amap.data, chunks=(1, 1))
    pair_map = sunpy.map.Map(da, amap.meta)
    assert isinstance(pair_map, sunpy.map.GenericMap)


def test_databaseentry():
    pytest.importorskip('sqlalchemy')
    sunpy_database = pytest.importorskip('sunpy.database')
    db = sunpy_database.Database(url='sqlite://', default_waveunit='angstrom')
    db.add_from_file(AIA_171_IMAGE)
    res = db.get_entry_by_id(1)
    db_map = sunpy.map.Map(res)
    assert isinstance(db_map, sunpy.map.GenericMap)


@pytest.mark.remote_data
def test_url_pattern():
    # A URL
    amap = sunpy.map.Map("http://data.sunpy.org/sample-data/AIA20110319_105400_0171.fits")
    assert isinstance(amap, sunpy.map.GenericMap)


def test_save():
    # Test save out
    aiamap = sunpy.map.Map(AIA_171_IMAGE)
    afilename = tempfile.NamedTemporaryFile(suffix='fits').name
    with pytest.warns(fits.verify.VerifyWarning, match="Invalid 'BLANK' keyword in header."):
        aiamap.save(afilename, filetype='fits', overwrite=True)
    backin = sunpy.map.Map(afilename)
    assert isinstance(backin, sunpy.map.sources.AIAMap)


@pytest.mark.remote_data
def test_map_list_urls_cache():
    """
    Test for https://github.com/sunpy/sunpy/issues/4006
    """
    urls = ['http://jsoc.stanford.edu/SUM80/D136597189/S00000/image_lev1.fits',
            'http://jsoc.stanford.edu/SUM79/D136597240/S00000/image_lev1.fits']
    sunpy.map.Map(urls)


@pytest.mark.filterwarnings('ignore:File may have been truncated')
@pytest.mark.parametrize('file, mapcls', [
    ["EIT_header/efz20040301.000010_s.header", sunpy.map.sources.EITMap],
    ["lasco_c2_25299383_s.header", sunpy.map.sources.LASCOMap],
    ["mdi.fd_Ic.20101015_230100_TAI.data.header", sunpy.map.sources.MDIMap],
    ["mdi.fd_M_96m_lev182.20101015_191200_TAI.data.header", sunpy.map.sources.MDIMap],
    ["euvi_20090615_000900_n4euA_s.header", sunpy.map.sources.EUVIMap],
    ["cor1_20090615_000500_s4c1A.header", sunpy.map.sources.CORMap],
    ["hi_20110910_114721_s7h2A.header", sunpy.map.sources.HIMap],
    [AIA_171_IMAGE, sunpy.map.sources.AIAMap],
    [RHESSI_IMAGE, sunpy.map.sources.RHESSIMap],
    ["FGMG4_20110214_030443.7.header", sunpy.map.sources.SOTMap],
    ["swap_lv1_20140606_000113.header", sunpy.map.sources.SWAPMap],
    ["HinodeXRT.header", sunpy.map.sources.XRTMap],
])
def test_sources(file, mapcls):
    p = pathlib.Path(get_test_filepath(file))
    if p.suffix == '.header':
        m = get_dummy_map_from_header(p)
    else:
        m = sunpy.map.Map(p)
    assert isinstance(m, mapcls)


def test_no_2d_hdus(tmpdir):
    # Create a fake FITS file with a valid header but 1D data
    tmp_fpath = str(tmpdir / 'data.fits')
    with fits.open(AIA_171_IMAGE, ignore_blank=True) as hdul:
        fits.writeto(tmp_fpath, np.arange(100), hdul[0].header)

    with pytest.raises(NoMapsInFileError, match='Found no HDUs with >= 2D data'):
        sunpy.map.Map(tmp_fpath)

    with pytest.warns(SunpyUserWarning, match='One of the arguments failed to parse'):
        sunpy.map.Map([tmp_fpath, AIA_171_IMAGE], silence_errors=True)


@skip_glymur
def test_map_jp2():
    jp2_map = sunpy.map.Map(AIA_193_JP2, memmap=False)
    assert isinstance(jp2_map, sunpy.map.GenericMap)
    # The base of an array that owns its memory is None
    # For some reason JP2 doesn't own its own memory but are not type of memmaps.
    assert jp2_map.data.base is not None
    jp2_map = sunpy.map.Map(AIA_193_JP2, memmap=True)
    assert isinstance(jp2_map, sunpy.map.GenericMap)
    assert jp2_map.data.base is not None


def test_map_fits():
    fits_map = sunpy.map.Map(AIA_171_IMAGE, memmap=False)
    assert isinstance(fits_map, sunpy.map.GenericMap)
    # The base of an array that owns its memory is None
    assert fits_map.data.base is None
    fits_map = sunpy.map.Map(AIA_171_IMAGE, memmap=True)
    assert isinstance(fits_map, sunpy.map.GenericMap)
    assert fits_map.data.base is not None
