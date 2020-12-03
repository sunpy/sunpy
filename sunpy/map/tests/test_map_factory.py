"""
Created on Fri Jun 21 15:05:09 2013

@author: stuart
"""
import os
import pathlib
import tempfile

import numpy as np
import pytest

from astropy.io import fits
from astropy.wcs import WCS

import sunpy
import sunpy.data.test
import sunpy.map
from sunpy.util.exceptions import NoMapsInFileError, SunpyUserWarning

filepath = pathlib.Path(sunpy.data.test.rootdir)
a_list_of_many = [os.fspath(f) for f in pathlib.Path(filepath, "EIT").glob("*")]
a_fname = a_list_of_many[0]


AIA_171_IMAGE = os.path.join(filepath, 'aia_171_level1.fits')
RHESSI_IMAGE = os.path.join(filepath, 'hsi_image_20101016_191218.fits')

amap = sunpy.map.Map(AIA_171_IMAGE)

valid_map_inputs = [(a_fname, ),
                    (pathlib.Path(a_fname), ),
                    (filepath / "EIT", ),
                    (os.fspath(filepath / "EIT"), ),
                    (filepath / "EIT" / "*", ),
                    (amap, ),
                    (amap.data, amap.meta),
                    ((amap.data, amap.meta), ),
                    ]


@pytest.mark.parametrize('args1', valid_map_inputs)
@pytest.mark.parametrize('args2', valid_map_inputs)
def test_two_map_inputs(args1, args2):
    out = sunpy.map.Map(*args1, *args2)
    if isinstance(out, list):
        for m in out:
            assert isinstance(m, sunpy.map.GenericMap)
    else:
        assert isinstance(out, sunpy.map.GenericMap)

# ==============================================================================
# Map Factory Tests
# ==============================================================================


class TestMap:
    def test_mapsequence(self):
        # Test making a MapSequence
        sequence = sunpy.map.Map(a_list_of_many, sequence=True)
        assert isinstance(sequence, sunpy.map.MapSequence)

    def test_composite(self):
        # Test making a CompositeMap
        comp = sunpy.map.Map(AIA_171_IMAGE, RHESSI_IMAGE, composite=True)
        assert isinstance(comp, sunpy.map.CompositeMap)

    # Want to check that patterns work, so ignore this warning that comes from
    # the AIA test data
    @pytest.mark.filterwarnings("ignore:Invalid 'BLANK' keyword in header")
    def test_patterns(self):
        # Test different Map pattern matching

        # File name
        eitmap = sunpy.map.Map(a_fname)
        assert isinstance(eitmap, sunpy.map.GenericMap)

        # Directory
        directory = pathlib.Path(filepath, "EIT")
        maps = sunpy.map.Map(os.fspath(directory))
        assert isinstance(maps, list)
        assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])
        # Test that returned maps are sorted
        files_sorted = sorted(list(directory.glob('*')))
        maps_sorted = [sunpy.map.Map(os.fspath(f)) for f in files_sorted]
        assert all([m.date == m_s.date for m, m_s in zip(maps, maps_sorted)])

        # Pathlib
        path = pathlib.Path(a_fname)
        eitmap = sunpy.map.Map(path)
        assert isinstance(eitmap, sunpy.map.GenericMap)
        maps = sunpy.map.Map(directory)
        assert isinstance(maps, list)
        assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])

        # Glob
        pattern = os.path.join(filepath, "EIT", "*")
        maps = sunpy.map.Map(pattern)
        assert isinstance(maps, list)
        assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])
        # Test that returned maps are sorted
        files_sorted = sorted(list(pathlib.Path(pattern).parent.glob('*')))
        maps_sorted = [sunpy.map.Map(os.fspath(f)) for f in files_sorted]
        assert all([m.date == m_s.date for m, m_s in zip(maps, maps_sorted)])
        # Single character wildcard (?)
        pattern = os.path.join(filepath, "EIT", "efz20040301.0?0010_s.fits")
        maps = sunpy.map.Map(pattern)
        assert isinstance(maps, list)
        assert len(maps) == 7
        assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])
        # Character ranges
        pattern = os.path.join(filepath, "EIT", "efz20040301.0[2-6]0010_s.fits")
        maps = sunpy.map.Map(pattern)
        assert isinstance(maps, list)
        assert len(maps) == 4
        assert ([isinstance(amap, sunpy.map.GenericMap) for amap in maps])

        # Already a Map
        amap = sunpy.map.Map(maps[0])
        assert isinstance(amap, sunpy.map.GenericMap)

        # A list of filenames
        maps = sunpy.map.Map(a_list_of_many)
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
        with fits.open(a_fname) as hdul:
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
        with pytest.warns(SunpyUserWarning, match='Missing CTYPE1 from metadata, assuming CTYPE1 is HPLN-TAN'):
            pair_map = sunpy.map.Map(data, header)
        assert isinstance(pair_map, sunpy.map.GenericMap)

        # Common keys not strings
        data = np.arange(0, 100).reshape(10, 10)
        header = {'cdelt1': 10, 'cdelt2': 10,
                  'telescop': 100,
                  'detector': 1,
                  'instrume': 50,
                  'cunit1': 'arcsec', 'cunit2': 'arcsec'}
        with pytest.warns(SunpyUserWarning, match='Missing CTYPE1 from metadata, assuming CTYPE1 is HPLN-TAN'):
            pair_map = sunpy.map.Map(data, header)
        assert isinstance(pair_map, sunpy.map.GenericMap)

    def test_errors(self, tmpdir):
        # If directory doesn't exist, make sure it's listed in the error msg
        nonexist_dir = 'nonexist'
        directory = pathlib.Path(filepath, nonexist_dir)
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

    # We want to check errors, so ignore warnings that are thrown
    @pytest.mark.filterwarnings("ignore:One of the data, header pairs failed to validate")
    @pytest.mark.parametrize('silence,error,match',
                             [(True, RuntimeError, 'No maps loaded'),
                              (False, sunpy.map.mapbase.MapMetaValidationError,
                               'Image coordinate units for axis 1 not present in metadata.')])
    def test_silence_errors(self, silence, error, match):
        # Check that the correct errors are raised depending on silence_errors value
        data = np.arange(0, 100).reshape(10, 10)
        header = {}
        with pytest.raises(error, match=match):
            pair_map = sunpy.map.Map(data, header, silence_errors=silence)

    # requires dask array to run properly
    def test_dask_array(self):
        dask_array = pytest.importorskip('dask.array')
        amap = sunpy.map.Map(AIA_171_IMAGE)
        da = dask_array.from_array(amap.data, chunks=(1, 1))
        pair_map = sunpy.map.Map(da, amap.meta)
        assert isinstance(pair_map, sunpy.map.GenericMap)

    # requires sqlalchemy to run properly
    def test_databaseentry(self):
        pytest.importorskip('sqlalchemy')
        sunpy_database = pytest.importorskip('sunpy.database')
        db = sunpy_database.Database(url='sqlite://', default_waveunit='angstrom')
        db.add_from_file(a_fname)
        res = db.get_entry_by_id(1)
        db_map = sunpy.map.Map(res)
        assert isinstance(db_map, sunpy.map.GenericMap)

    @pytest.mark.remote_data
    def test_url_pattern(self):
        # A URL
        amap = sunpy.map.Map("http://data.sunpy.org/sample-data/AIA20110319_105400_0171.fits")
        assert isinstance(amap, sunpy.map.GenericMap)

    def test_save(self):
        # Test save out
        eitmap = sunpy.map.Map(a_fname)
        afilename = tempfile.NamedTemporaryFile(suffix='fits').name
        with pytest.warns(SunpyUserWarning, match='The meta key  is not valid ascii'):
            eitmap.save(afilename, filetype='fits', overwrite=True)
        backin = sunpy.map.Map(afilename)
        assert isinstance(backin, sunpy.map.sources.EITMap)


@pytest.mark.remote_data
def test_map_list_urls_cache():
    """
    Test for https://github.com/sunpy/sunpy/issues/4006
    """
    urls = ['http://jsoc.stanford.edu/SUM80/D136597189/S00000/image_lev1.fits',
            'http://jsoc.stanford.edu/SUM79/D136597240/S00000/image_lev1.fits']

    sunpy.map.Map(urls)


# TODO: Test HMIMap, SXTMap
#
# Catch Hinode/XRT warning
@pytest.mark.filterwarnings('ignore:File may have been truncated')
@pytest.mark.parametrize('file, mapcls',
                         [[filepath / 'EIT' / "efz20040301.000010_s.fits", sunpy.map.sources.EITMap],
                          [filepath / "lasco_c2_25299383_s.fts", sunpy.map.sources.LASCOMap],
                          [filepath / "mdi_fd_Ic_6h_01d.5871.0000_s.fits", sunpy.map.sources.MDIMap],
                          [filepath / "mdi_fd_M_96m_01d.5874.0005_s.fits", sunpy.map.sources.MDIMap],
                          [filepath / "euvi_20090615_000900_n4euA_s.fts", sunpy.map.sources.EUVIMap],
                          [filepath / "cor1_20090615_000500_s4c1A.fts", sunpy.map.sources.CORMap],
                          [filepath / "hi_20110910_114721_s7h2A.fts", sunpy.map.sources.HIMap],
                          [AIA_171_IMAGE, sunpy.map.sources.AIAMap],
                          [RHESSI_IMAGE, sunpy.map.sources.RHESSIMap],
                          [filepath / "FGMG4_20110214_030443.7.fits", sunpy.map.sources.SOTMap],
                          [filepath / "swap_lv1_20140606_000113.fits", sunpy.map.sources.SWAPMap],
                          [filepath / "HinodeXRT.fits", sunpy.map.sources.XRTMap]
                          ]
                         )
def test_sources(file, mapcls):
    m = sunpy.map.Map(file)
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
