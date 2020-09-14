"""Test cases for SOHO Map subclasses.
This particular test file pertains to MDIMap.
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
from sunpy.coordinates import frames
from sunpy.map import Map
from sunpy.map.sources.soho import MDIMap, MDISynopticMap
from sunpy.util.exceptions import SunpyMetadataWarning


@pytest.fixture
def mdi():
    path = sunpy.data.test.rootdir
    fitspath = glob.glob(os.path.join(path, "mdi_fd_Ic_6h_01d.5871.0000_s.fits"))
    return Map(fitspath)


@pytest.fixture
def mdi_synoptic():
    # sample MDI map
    # Please do not edit this header, it is a 1:1 copy of a MDI map header
    header = dedent("""
        SIMPLE  =                    T / file does conform to FITS standard
        BITPIX  =                  -32 / number of bits per data pixel
        NAXIS   =                    2 / number of data axes
        NAXIS1  =                  720 / length of data axis 1
        NAXIS2  =                  360 / length of data axis 2
        EXTEND  =                    T / FITS dataset may contain extensions
        COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
        COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
        DATE    = '2017-05-17T00:00:00'
        TELESCOP= 'SOHO'
        INSTRUME= 'MDI'
        WAVELNTH= 6768.000000
        BUNIT   = 'Mx/cm^2'
        CONTENT = 'Carrington Synoptic Chart Of Mr Field'
        HISTORY Carrington-Time conversion corrected; o2helio.c bug corrected -- July 20
        HISTORY 13
        BLD_VERS=          -2147483648
        CTYPE1  = 'CRLN-CEA'
        CTYPE2  = 'CRLT-CEA'
        CRPIX1  = 360.40
        CRPIX2  = 180.50
        CRVAL1  = 755460.00
        CRVAL2  = 0.00
        CDELT1  = -0.500000
        CDELT2  = 0.006944
        CUNIT1  = 'Degree'
        CUNIT2  = 'Sine Latitude'
        CRDER1  = 'nan     '
        CRDER2  = 'nan     '
        CSYSER1 = 'nan     '
        CSYSER2 = 'nan     '
        WCSNAME = 'Carrington Heliographic'
        T_REC   = '-4712.01.01_11:59:28_TAI'
        TRECEPOC= '1977.01.01_00:00:00_TAI'
        TRECSTEP= 2356586.000000
        TRECUNIT= 'secs'
        CADENCE = 360.000000
        DATASIGN=                    1
        T_OBS   = '2010.07.27_00:09:05_TAI'
        T_START = '2010.07.13_09:38:38_TAI'
        T_STOP  = '2010.08.09_14:48:26_TAI'
        T_ROT   = '2010.07.27_00:09:05_TAI'
        CAR_ROT =                 2099
        CARRTIME= 755460.00
        B0_ROT  = 5.379462
        B0_FRST = 4.150289
        B0_LAST = 6.337214
        LON_FRST= 755280.20
        LON_LAST= 755639.70
        LON_STEP= -0.500000
        W_OFFSET= 0.000000
        W_WEIGHT= 'Even'
        IMG_NUM =                  338
        IMG_FRST= '2010.07.12_00:00:00_TAI'
        IMG_LAST= '2010.08.10_00:00:00_TAI'
        IMG_ROT = '2010.07.24_06:24:00_TAI'
        HWNWIDTH= 20.000000
        EQPOINTS= 20.000000
        NSIGMA  = 3.000000
        CARSTRCH=                    1
        DIFROT_A= 13.562000
        DIFROT_B= -2.040000
        DIFROT_C= -1.487500
        TOTVALS =               259200
        DATAVALS=               257588
        MISSVALS=                 1612
        DATAMIN = -1880.387695
        DATAMAX = 1857.843750
        DATAMEDN= 0.072514
        DATAMEAN= 0.025515
        DATARMS = 26.509175
        DATASKEW= -10.533994
        DATAKURT= 1172.076612
        CALVER64=                    0
        RECNUM  =                  188
        CHECKSUM= '8MHfAK9d2KGd8K9d'   / HDU checksum updated 2020-04-23T08:48:07
        DATASUM = '3044927737'         / data unit checksum updated 2020-04-23T08:48:07""")
    header = fits.Header.fromstring(header, sep='\n')
    data = np.random.rand(header['naxis1'], header['naxis2'])
    return Map(data, header)


# MDI Tests
def test_fitstoMDI(mdi):
    """Tests the creation of MDIMap using FITS."""
    assert isinstance(mdi, MDIMap)


def test_is_datasource_for(mdi):
    """Test the is_datasource_for method of MDIMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert mdi.is_datasource_for(mdi.data, mdi.meta)


def test_observatory(mdi):
    """Tests the observatory property of the MDIMap object."""
    assert mdi.observatory == "SOHO"


def test_measurement(mdi):
    """Tests the measurement property of the MDIMap object."""
    assert mdi.measurement == "continuum"


def test_waveunit(mdi):
    assert mdi.waveunit == "Angstrom"


def test_observer(mdi):
    assert isinstance(mdi.observer_coordinate.frame, frames.HeliographicStonyhurst)
    assert u.allclose(mdi.observer_coordinate.lat, -5.774028172878*u.deg)
    assert u.allclose(mdi.observer_coordinate.lon, -0.10522355*u.deg)
    assert u.allclose(mdi.observer_coordinate.radius, 0.9739569156244*u.AU)


def test_carrington(mdi):
    assert u.allclose(mdi.carrington_longitude, mdi.meta['obs_l0']*u.deg)
    assert u.allclose(mdi.carrington_latitude, mdi.meta['obs_b0']*u.deg)


@pytest.mark.filterwarnings("error")
def test_synoptic_source(mdi_synoptic):
    assert isinstance(mdi_synoptic, MDISynopticMap)
    # Check that the WCS is valid
    with pytest.warns(SunpyMetadataWarning, match='Missing metadata for observer'):
        mdi_synoptic.wcs
