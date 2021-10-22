from textwrap import dedent

import numpy as np
import pytest

import astropy.units as u
from astropy.io import fits

from sunpy.map import Map
from sunpy.map.sources.hinode import SOTMap
from sunpy.util.exceptions import SunpyMetadataWarning


@pytest.fixture
def sot():
    raw_header = dedent("""\
        SIMPLE  =                    T / created Wed Oct 21 17:17:46 2015
        BITPIX  =                   16
        NAXIS   =                    2
        NAXIS1  =                 2048
        NAXIS2  =                 1024
        EXTEND  =                    T
        DATE    = '2015-10-21T17:17:46.000'        /creation date
        DATE_RF0= '2015-10-21T17:17:46.000'        /creation date
        TELESCOP= 'HINODE'
        MDP_CLK =            361702812
        FILEORIG= '2015_1021_171239.sci'
        MDPCTREF=            361701639
        CTREF   =           1801875595
        CTRATE  =              584.000
        TIMEERR =                    0
        EXP0    =             0.123088
        OBT_TIME=            361701866
        OBT_END =            361701929
        DATE_OBS= '2015-10-13T23:13:44.601'
        TIME-OBS= '23:13:44.601'
        CTIME   = 'Tue Oct 13 23:13:44 2015'
        DATE_END= '2015-10-13T23:13:44.724'
        TIMESPAN=             0.123088
        TIMESYS = 'UTC'
        INSTRUME= 'SOT/WB'
        ORIGIN  = 'JAXA/ISAS, SIRIUS'
        DATA_LEV=                    0
        ORIG_RF0= 'JAXA/ISAS, SIRIUS'
        VER_RF0 = '1.66'
        PROG_VER=                  609
        SEQN_VER=                  919
        PARM_VER=                  628
        PROG_NO =                   17
        SUBR_NO =                    1
        SEQN_NO =                   27
        MAIN_CNT=                    1
        MAIN_RPT=                    1
        MAIN_POS=                    1
        SUBR_CNT=                    1
        SUBR_RPT=                    2
        SUBR_POS=                    1
        SEQN_CNT=                    1
        SEQN_RPT=                    1
        SEQN_POS=                    5
        OBSTITLE= ' '
        TARGET  = ' '
        SCI_OBJ = ' '
        SCI_OBS = 'TBD'
        OBS_DEC = ' '
        JOIN_SB = ' '
        OBS_NUM =                    0
        JOP_ID  =                    0
        NOAA_NUM=                    0
        OBSERVER= '  '
        PLANNER = '  '
        TOHBANS = '  '
        DATATYPE= 'SCI'
        FLFLG   = 'NON'
        OBS_MODE= 'QT'
        SAA     = 'OUT'
        HLZ     = 'OUT'
        OBS_ID  =                    1
        GEN_ID  =                    1
        FRM_ID  =                    2
        WAVEID  =                    2
        OBS_TYPE= 'FG (simple)'
        MACROID =                13348
        XSCALE  =             0.108960
        YSCALE  =             0.108960
        FGXOFF  =                    0
        FGYOFF  =                    0
        FGCCDIX0=                    0
        FGCCDIX1=                 4095
        FGCCDIY0=                    0
        FGCCDIY1=                 2047
        CRPIX1  =              1024.50
        CRPIX2  =              512.500
        SC_ATTX =             -15.8358
        SC_ATTY =              19.2347
        CRVAL1  =             -15.8358
        CRVAL2  =              19.2347
        CDELT1  =             0.108960
        CDELT2  =             0.108960
        CUNIT1  = 'arcsec'
        CUNIT2  = 'arcsec'
        CTYPE1  = 'Solar-X'
        CTYPE2  = 'Solar-Y'
        SAT_ROT =              0.00000
        INST_ROT=             0.412000
        CROTA1  =             0.412000
        CROTA2  =             0.412000
        XCEN    =             -15.8358
        YCEN    =              19.2347
        FOVX    =              223.150
        FOVY    =              111.575
        TR_MODE = 'FIX'
        FGBINX  =                    1
        FGBINY  =                    1
        EXPTIME =             0.122880
        WAVE    = 'Ca II H line'
        DARKFLAG=                    0
        BITCOMP1=                    6
        IMGCOMP1=                    7
        QTABLE1 =                    2
        BITCOMP2=                    1
        IMGCOMP2=                    7
        QTABLE2 =                    1
        PCK_SN0 =             42477799
        PCK_SN1 =             42477831
        NUM_PCKS=                   33
        FGMODE  = 'shuttered'
        FGNINT  =                    1
        ROILOOP =                    0
        NROILOOP=                    0
        CTSERVO =                    1
        CTMESTAT=                36864
        CTMEX   =                 -801
        CTMEY   =                  329
        CTMODE  =                   35
        T_SPCCD =             -43.4724
        T_FGCCD =             -32.2096
        T_CTCCD =              26.6964
        T_SPCEB =             -6.63403
        T_FGCEB =              1.65425
        T_CTCEB =              2.92207
        MASK    =                   22
        WBFW    =                  118
        WEDGE   =                   21
        NBFW    =                  118
        TF1     =                  132
        TF2     =                  133
        TF3     =                  137
        TF4     =                  129
        TF5     =                  106
        TF6     =                   54
        TF7     =                   17
        TF8     =                  154
        SLITENC =                 1991
        FOCUS   =                 2035
        WBEXP   =                   50
        NBEXP   =                  199
        WAVEOFF =                    0
        ROISTART=                    0
        ROISTOP =                 1025
        DOPVUSED=                 -110
        CAMGAIN =                    2
        CAMDACA =                    8
        CAMDACB =                    8
        CAMPSUM =                    2
        CAMSSUM =                    2
        CAMAMP  =                    0
        CAMSCLK =                    0
        PMUDELAY=                  128
        BITCVER1=                45094
        ACHFVER1=                40961
        DCHFVER1=                53249
        QTABVER1=                57365
        BYTECNTI=               401402
        PIXCNTI =              2097152
        BITSPPI =              1.53123
        BYTECNTQ=                    0
        PIXCNTQ =                    0
        BITSPPQ =              0.00000
        PERCENTD=              100.000
        """)
    header = fits.Header.fromstring(raw_header, sep='\n')
    data = np.random.rand(header['naxis1'], header['naxis2'])
    return Map(data, header)


# SOT Tests
def test_fitstoSOT(sot):
    """Tests the creation of SOTMap using FITS."""
    assert isinstance(sot, SOTMap)


def test_is_datasource_for(sot):
    """Test the is_datasource_for method of SOTMap.
    Note that header data to be provided as an argument
    can be a MetaDict object."""
    assert sot.is_datasource_for(sot.data, sot.meta)


def test_observatory(sot):
    """Tests the observatory property of the SOTMap object."""
    assert sot.observatory == "Hinode"


def test_measurement(sot):
    """Tests the measurement property of the SOTMap object."""
    assert sot.measurement is None


def test_instruments(sot):
    """Tests the Instruments object of SOTMap."""
    assert (sot.Instruments == ['SOT/WB',
                                'SOT/NB', 'SOT/SP', 'SOT/CT'])


def test_waves(sot):
    """Tests the Waves object of SOTMap."""
    assert (sot.Waves == ['6302A', 'BFI no move',
                          'CN bandhead 3883', 'Ca II H line',
                          'G band 4305', 'NFI no move', 'TF Fe I 6302',
                          'TF Mg I 5172', 'TF Na I 5896',
                          'blue cont 4504', 'green cont 5550',
                          'red cont 6684'])


def test_obstype(sot):
    """Tests the Observation_Type object of SOTMap."""
    assert (sot.Observation_Type == ['FG (simple)',
                                     'FG focus scan', 'FG shuttered I and V',
                                     'FG shutterless I and V', 'FG shutterless I and V with 0.2s intervals',
                                     'FG shutterless Stokes', 'SP IQUV 4D array'])


def test_wcs(sot):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    with pytest.warns(SunpyMetadataWarning, match='assuming Earth-based observer'):
        sot.pixel_to_world(0*u.pix, 0*u.pix)
