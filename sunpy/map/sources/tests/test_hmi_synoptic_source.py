from textwrap import dedent

import numpy as np
import pytest

from astropy.io import fits

from sunpy.map import Map
from sunpy.map.sources.sdo import HMISynopticMap


@pytest.fixture
def hmi_synoptic():
    # sample HMI synoptic map
    # Please do not edit this header, it is a 1:1 copy of a HMI synoptic map header
    raw_header = dedent("""\
        SIMPLE  =                    T / file does conform to FITS standard
        BITPIX  =                  -32 / number of bits per data pixel
        NAXIS   =                    2 / number of data axes
        NAXIS1  =                  720 / length of data axis 1
        NAXIS2  =                  360 / length of data axis 2
        EXTEND  =                    T / FITS dataset may contain extensions
        COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
        COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H
        DATE    = '2018-11-29T00:00:00'
        TELESCOP= 'SDO/HMI'
        INSTRUME= 'HMI_SIDE1'
        WAVELNTH= 6173.000000
        BUNIT   = 'Mx/cm^2'
        CONTENT = 'Carrington Synoptic Chart Of Br Field'
        HISTORY Carrington-Time conversion corrected; o2helio.c bug corrected -- July 20
        HISTORY 13
        BLD_VERS=          -2147483648
        CTYPE1  = 'CRLN-CEA'
        CTYPE2  = 'CRLT-CEA'
        CRPIX1  = 360.400000
        CRPIX2  = 180.500000
        CRVAL1  = 795420.000000
        CRVAL2  = 0.000000
        CDELT1  = -0.500000
        CDELT2  = 0.005556
        CUNIT1  = 'Degree'
        CUNIT2  = 'Sine Latitude'
        CRDER1  = 'nan     '
        CRDER2  = 'nan     '
        CSYSER1 = 'nan     '
        CSYSER2 = 'nan     '
        WCSNAME = 'Carrington Heliographic'
        CADENCE = 360.000000
        DATASIGN=                    1
        T_OBS   = '2018.11.09_12:30:52_TAI'
        T_START = '2018.10.26_20:53:36_TAI'
        T_STOP  = '2018.11.23_04:13:10_TAI'
        T_ROT   = '2018.11.09_12:30:52_TAI'
        T_EARTH = ''
        CAR_ROT =                 2210
        CARRTIME= 795420.000000
        B0_ROT  = 3.504195
        B0_FRST = 4.903270
        B0_LAST = 1.897488
        EARTH_B0= 'nan     '
        LON_FRST= 795240.200000
        LON_LAST= 795599.700000
        LON_STEP= -0.500000
        W_OFFSET= 0.000000
        W_WEIGHT= 'Even'
        IMG_NUM =                 3276
        IMG_FRST= '2018.10.26_00:00:00_TAI'
        IMG_LAST= '2018.11.23_23:00:00_TAI'
        IMG_ROT = '2018.11.09_12:36:00_TAI'
        HWNWIDTH= 15.000000
        EQPOINTS= 20.000000
        NSIGMA  = 3.000000
        CARSTRCH=                    1
        DIFROT_A= 13.562000
        DIFROT_B= -2.040000
        DIFROT_C= -1.487500
        TOTVALS =               259200
        DATAVALS=               258152
        MISSVALS=                 1048
        DATAMIN = -275.919403
        DATAMAX = 352.637482
        DATAMEDN= 0.006980
        DATAMEAN= 0.024495
        DATARMS = 8.102401
        DATASKEW= 0.331222
        DATAKURT= 70.596480
        CALVER64=               270610
        CODEVER = ''
        QUALITY =                    0
        RECNUM  =                  232
        CHECKSUM= 'cjHPdiENciENciEN'   / HDU checksum updated 2020-04-22T21:01:56
        DATASUM = '2978664836'         / data unit checksum updated 2020-04-22T21:01:56
        """)

    header = fits.Header.fromstring(raw_header, sep='\n')
    data = np.random.rand(header['naxis1'], header['naxis2'])
    return Map(data, header)


def test_fitstoHMISynoptic(hmi_synoptic):
    """Tests the creation of HMISynopticMap using FITS."""
    assert isinstance(hmi_synoptic, HMISynopticMap)


def test_is_datasource_for(hmi_synoptic):
    """Test the is_datasource_for method of HMISynopticMap.
    Note that header data to be provided as an argument
    can be a MetaDict object, which in this case is
    hmi.meta."""
    assert hmi_synoptic.is_datasource_for(hmi_synoptic.data, hmi_synoptic.meta)


def test_observatory(hmi_synoptic):
    """Tests the observatory property of the HMISynopticMap object."""
    assert hmi_synoptic.observatory == "SDO"


def test_measurement(hmi_synoptic):
    """Tests the measurement property of the HMISynopticMap object."""
    assert hmi_synoptic.measurement == "carrington"
