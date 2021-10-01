"""Test cases for SOHO Map subclasses.
This particular test file pertains to LASCOMap.
@Author: Pritish C. (VaticanCameos)
"""

import os
import glob
from textwrap import dedent

import numpy as np
import pytest

import astropy.units as u
from astropy.io import fits

import sunpy.data.test
from sunpy.map import Map
from sunpy.map.sources.soho import LASCOMap
from sunpy.tests.helpers import skip_glymur
from sunpy.time import parse_time
from sunpy.util.exceptions import SunpyMetadataWarning

path = sunpy.data.test.rootdir


def lascoc3():
    # Please do not edit this header, it is a 1:1 copy of a LASCO C3 header
    raw_header = dedent("""\
        SIMPLE  =                    T / Written by IDL:  Thu Jun  6 19:03:59 2002
        BITPIX  =                   16 /
        NAXIS   =                    2 /
        NAXIS1  =                 1024 /
        NAXIS2  =                 1024 /
        FILENAME= '32088304.fts'       /
        FILEORIG= '020521_001853.img'  /
        DATE    = '2002/06/06 23:03:55.204' /
        DATE-OBS= '2002/05/21'         /
        TIME-OBS= '00:18:06.516'       /
        P1COL   =                   20 /
        P1ROW   =                    1 /
        P2COL   =                 1043 /
        P2ROW   =                 1024 /
        VERSION =                    2 /
        EXPTIME =              19.0996 /
        EXP0    =              17.0000 /
        EXPCMD  =              17.0000 /
        EXP1    =              1.88135 /
        EXP2    =              3.17676 /
        EXP3    =              2.09961 /
        TELESCOP= 'SOHO    '           /
        INSTRUME= 'LASCO   '           /
        DETECTOR= 'C3      '           /
        READPORT= 'C       '           /
        SUMROW  =                    0 /
        SUMCOL  =                    0 /
        LEBXSUM =                    1 /
        LEBYSUM =                    1 /
        SHUTTR  =                    0 /
        LAMP    =                    0 /
        FILTER  = 'Clear   '           /
        POLAR   = 'Clear   '           /
        LP_NUM  = 'Normal  '           /
        OS_NUM  =                 3390 /
        IMGCTR  =                  717 /
        IMGSEQ  =                    0 /
        COMPRSSN= 'X#      '           /
        HCOMP_SF=                   32 /
        MID_DATE=                52415 /
        MID_TIME=              1096.07 /
        PLATESCL=              56.0000 /
        OFFSET  =              378.876 /
        IMAGE_CT=                  717 /
        SEQ_NUM =                    0 /
        OBT_TIME=        1.4006315E+09 /
        R1COL   =                   20 /
        R1ROW   =                    1 /
        R2COL   =                 1043 /
        R2ROW   =                 1024 /
        EFFPORT = 'C       '           /
        RECTIFY = 'TRUE    '           /
        DATAMIN =              376.000 /
        DATAMAX =              11183.0 /
        DATAZER =               101375 /
        DATASAT =                    0 /
        DATAAVG =              1886.68 /
        DATASIG =              1630.27 /
        DATAP01 =                  396 /
        DATAP10 =                  780 /
        DATAP25 =                  968 /
        DATAP75 =                 2029 /
        DATAP90 =                 4213 /
        DATAP95 =                 5959 /
        DATAP98 =                 7190 /
        DATAP99 =                 8023 /
        CRPIX1  =            517.95599 /
        CRPIX2  =            532.63202 /
        CRVAL1  =              0.00000 /
        CRVAL2  =              0.00000 /
        CROTA1  =              0.00000 /
        CROTA2  =              0.00000 /
        CTYPE1  = 'SOLAR-X '           /
        CTYPE2  = 'SOLAR-Y '           /
        CUNIT1  = 'ARCSEC  '           /
        CUNIT2  = 'ARCSEC  '           /
        CDELT1  =            56.000000 /
        CDELT2  =            56.000000 /
        HISTORY offset_bias.pro	1.24 12/13/01, 378.876
        HISTORY V15 18 Jun 1997 MAKE_FITS_HDR
        HISTORY 2
        """)
    header = fits.Header.fromstring(raw_header, sep='\n')
    data = np.random.rand(header['naxis1'], header['naxis2'])
    return Map(data, header)


def lasco():
    fitspath = glob.glob(os.path.join(path, "lasco_c2_25299383_s.fts"))
    return Map(fitspath)


@pytest.fixture(params=[lasco, lascoc3], ids=['C2', 'C3'])
def lasco_map(request):
    return request.param()


@pytest.fixture
def lasco_helioviewer():
    jp2path = glob.glob(os.path.join(
        path, "2013_05_13__16_54_06_137__SOHO_LASCO_C3_white-light.jp2"))
    return Map(jp2path)


def test_fitstoLASCO(lasco_map):
    """Tests the creation of LASCOMap using FITS."""
    assert isinstance(lasco_map, LASCOMap)


def test_is_datasource_for(lasco_map):
    """Test the is_datasource_for method of LASCOMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert lasco_map.is_datasource_for(lasco_map.data, lasco_map.meta)


def test_measurement(lasco_map):
    """Tests the measurement property of the LASCOMap object."""
    assert lasco_map.measurement == "white-light"


def test_wavelength(lasco_map):
    """Tests wavelength property."""
    assert lasco_map.wavelength is None


def test_date(lasco_map):
    assert lasco_map.date == parse_time(
        {'C2': '2009-02-28T00:05:33.380',
         'C3': '2002-05-21T00:18:06.516'}[lasco_map.detector])


def test_nickname(lasco_map):
    assert lasco_map.nickname == {'C2': 'LASCO-C2 Orange',
                                  'C3': 'LASCO-C3 Clear'}[lasco_map.detector]


def test_observatory(lasco_map):
    """Tests the observatory property of the LASCOMap object."""
    assert lasco_map.observatory == "SOHO"


def test_norm_clip(lasco_map):
    # Tests that the default normalizer has clipping disabled
    assert not lasco_map.plot_settings['norm'].clip


@skip_glymur
def test_helioviewer_rotation(lasco_map, lasco_helioviewer):
    """Tests that rotation metadata is correctly removed
    for JPEG2000 images provided by Helioviewer.org."""
    rmatrix = {'C2': [[0.999966, -0.008296], [0.008296, 0.999966]],
               'C3': [[1, 0], [0, 1]]}[lasco_map.detector]
    np.testing.assert_allclose(lasco_map.rotation_matrix, rmatrix, rtol=1e-6)
    np.testing.assert_array_equal(lasco_helioviewer.rotation_matrix, [[1., 0.], [0., 1.]])


def test_wcs(lasco_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        lasco_map.pixel_to_world(0*u.pix, 0*u.pix)
