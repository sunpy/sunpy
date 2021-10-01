"""Tests for EUI Solar Orbiter Map"""
from textwrap import dedent

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

from sunpy.map import Map
from sunpy.map.sources import EUIMap


@pytest.fixture
def eui_fsi_map(scope="module"):
    # Please do not edit this header, it is a 1:1 copy of a SOLO EUI FSI 304 map
    raw_header = dedent("""\
        SIMPLE  =                    T / conforms to FITS standard
        BITPIX  =                   16 / array data type
        NAXIS   =                    2 / number of array dimensions
        NAXIS1  =                  768
        NAXIS2  =                  768
        LONGSTRN= 'OGIP 1.0'           / The OGIP long string convention may be used
        FILENAME= 'solo_L1_eui-fsi304-image_20201021T145510206_V01.fits' / FITS filename
        DATE    = '2020-11-23T16:40:54.714' / [UTC] FITS file creation date
        FILE_RAW= 'BatchRequest.PktTmRaw.SOL.0.2020.295.15.15.01.857.eJeU@2020.295.15.&'
        CONTINUE  '15.03.463.1.xml&'
        CONTINUE  '' / raw filename
        PARENT  = 'solo_L0_eui-fsi###-image_0656607273e84f_V00.fits' / source file curre
        APID    =                  924 / APID number of associated TM
        DATE-OBS= '2020-10-21T14:55:10.206' / [UTC] deprecated, same as DATE-BEG.
        DATE-BEG= '2020-10-21T14:55:10.206' / [UTC] start time of observation
        DATE-AVG= '2020-10-21T14:55:13.206' / [UTC] average time of observation
        TIMESYS = 'UTC     '           / system used for time keywords
        OBT_BEG =    656607273.9074554 / start acquisition time in OBT
        LEVEL   = 'L1      '           / data processing level
        ORIGIN  = 'Royal Observatory of Belgium' / file generation location
        VERS_SW = '807     '           / version of SW that provided FITS file
        VERSION = '01      '           / incremental version number
        IMGTYPE = 'solar image'        / type of image; solar, calib., engineering
        JOBID   = '20201123T163957.998Z_01f' / unique pipeline job ID
        COMPLETE= 'C       '           / C or I for complete/incomplete
        OBSRVTRY= 'Solar Orbiter'      / satellite name
        TELESCOP= 'SOLO/EUI/FSI'       / telescope/Sensor name
        INSTRUME= 'EUI     '           / instrument name
        DETECTOR= 'FSI     '           / instrument subunit or sensor
        DOORINT = 'open    '           / internal door position
        DOOREXT = 'open    '           / external HS Door 2 FSI @-76.08 s
        XPOSURE =                  6.0 / [s] total effective exposure time
        FILTER  = 'Magnesium_304_2'    / filter position
        WAVELNTH=                  304 / [Angstrom] characteristic wavelength observatio
        WAVEMIN =                  250 / [Angstrom] min wavelength resp. > 0.05 of max
        WAVEMAX =                  350 / [Angstrom] max wavelength resp. > 0.05 of max
        SOOPNAME= 'not defined'        / name of the SOOP Campaign that the data belong
        SOOPTYPE= '000     '           / campaign ID(s) that the data belong to
        OBS_MODE= 'GENERIC_CALIB'      / observation mode
        OBS_TYPE= '2ZpG    '           / encoded version of OBS_MODE
        OBS_ID  = 'SEUI_021A_000_000_2ZpG_11K' / unique ID of the individual observation
        TARGET  = 'not defined'        / type of target from planning
        BTYPE   = 'Intensity'          / type of data
        BUNIT   = 'DN      '           / units of physical value, after BSCALE, BZERO
        UCD     = 'phot.count;em.UV.10-50nm' / Unified Content Descriptor
        BLANK   =                32767 / value undefined pixels before BSCALE,BZERO
        PXBEG1  =                    1 / first read-out pixel in dimension 1
        PXEND1  =                 3072 / last read-out pixel in dimension 1
        PXBEG2  =                    1 / first read-out pixel in dimension 2
        PXEND2  =                 3072 / last read-out pixel in dimension 2
        NBIN1   =                    4 / on-board binning factor in dim 1
        NBIN2   =                    4 / on-board binning factor in dim 1
        NBIN    =                   16 / total binning factor
        WCSNAME = 'Helioprojective-cartesian' / name of coordinate system
        CTYPE1  = 'HPLN-TAN'           / helioprojective longitude (Solar X)
        CTYPE2  = 'HPLT-TAN'           / helioprojective latitude (Solar Y)
        CUNIT1  = 'arcsec  '           / units along axis 1
        CUNIT2  = 'arcsec  '           / units along axis 2
        PC1_1   =   0.9999102177627085 / WCS coordinate transformation matrix
        PC1_2   = -0.01342635032202125 / WCS coordinate transformation matrix
        PC2_1   =  0.01337343428604044 / WCS coordinate transformation matrix
        PC2_2   =   0.9999102177627085 / WCS coordinate transformation matrix
        CDELT1  =             17.74456 / [arcsec] pixel scale along axis 1
        CDELT2  =           17.7796312 / [arcsec] pixel scale along axis 2
        CROTA   =   0.7677787561176675 / [deg] rotation angle
        CRVAL1  =    111.3320447860757 / [arcsec] value of reference pixel along axis 1
        CRVAL2  =    112.4105857957075 / [arcsec] value of reference pixel along axis 2
        CRPIX1  =                384.5 / [pixel] reference pixel location along axis 1
        CRPIX2  =                384.5 / [pixel] reference pixel location along axis 2
        LONPOLE =                180.0 / [deg] native longitude of the celestial pole
        ATT_SKD = 'v106_20201120_002+flown' / attitude SKD version, +flown/+predicted
        DETGAINL=                0.875 / commanded low gain value
        DETGAINH=                  3.0 / commanded high-gain value
        GAINCOMB= 'combined'           / commanded low/high gain combination
        READOUTM=                    4 / commanded FEE readout mode
        DOWNLOAM=                    0 / commanded FEE download mode
        GAINTHRE=                15600 / commanded threshold value for H/L gain
        LEDSTATE= 'all off '           / LED control for current telescope
        TEMPINT =     211.159425070688 / [K] internal APS detector temperature
        TEMP1DET=    211.3479264844486 / [K] last measured APS det. T before date-avg
        TEMP2DET=     211.159425070688 / [K] earliest measured APS det. T after date-avg
        TTEMP1  = '2020-10-21T14:55:00.902583' / [UTC] time TEMP1DET measurement
        TTEMP2  = '2020-10-21T14:56:00.901291' / [UTC] time TEMP2DET measurement
        ALU     =                21131 / CEB ALU register
        ALU2    =                    1 / CEB ALU2 register
        DETREG20=                   85 / REG20_VDAC_CLIPP
        DETREG21=                   87 / REG21_VDAC_OFFSETP
        DETREG22=                   63 / REG22_VDAC_CMREF
        DETREG23=                   39 / REG23_VDAC_OFFSETN
        DETREG24=                   41 / REG24_VDAC_CLIPN
        DETREG25=                   35 / REG25_VDAC_CMREF_LV
        DETREG26=                  137 / REG26_IDAC_CDSSTAGE2_3
        DETREG27=                  136 / REG27_IDAC_CDSSTAGE1_COMPA
        DETREG28=                  136 / REG28_IDAC_INVBUFFER_REFBU
        DETREG29=                  136 / REG29_IDAC_COLBUFFER_COLBU
        DETREG2A=                    8 / REG2A_IDAC_COLPC_COLGAINST
        DETREG2B=                  136 / REG2B_IDAC_OUTPUTDRIVER_CO
        DETREG2C=                    0 / REG2C_VDAC_BLACKSUN_EVEN
        DETREG2D=                  136 / REG2D_IDAC_ABS_REBUFF_TEMP
        DETREG2E=                  255 / REG2E_COLGAIN_EVEN_FF_MID
        DETREG2F=                  240 / REG2F_COLGAIN_EVEN_FF_LOW
        DETREG30=                    0 / REG30_COLGAIN_EVEN_FB_MID
        DETREG31=                   15 / REG31_COLGAIN_EVEN_FB_LOW
        DETREG32=                    0 / REG32_COLGAIN_EVEN_FB_HIGH
        DETREG33=                  127 / REG33_COLGAIN_ODD_FF_MID
        DETREG34=                    0 / REG34_COLGAIN_ODD_FF_LOW
        DETREG35=                    0 / REG35_COLGAIN_ODD_FB_MID
        DETREG36=                  255 / REG36_COLGAIN_ODD_FB_LOW
        DETREG37=                    0 / REG37_COLGAIN_VDAC_SIGCLAM
        DETREG38=                   61 / REG38_CDS_EN_SAMPLE_CLOCK_
        DETREG39=                    0 / REG39_MBS_PIXCOL_ADDR_LOW
        DETREG3A=                    0 / REG3A_MBS_PIXCOL_ADDR_HIGH
        DETREG3B=                    0 / REG3B_MBS_MUXBUS_SR_EOSX_S
        DETREG3C=                    0 / REG3C_VDAC_SIGC_LAMP_BLACK
        DETREG3D=                    0 / REG3D_XWIN_ADDRESS
        DETREG3E=                   65 / REG3E_VDAC_BUSCLAMPHIGH
        DETREG3F=                   65 / REG3F_VDAC_BUSCLAMPLOW
        DOORPOS =                   34 / Door position (raw) = open
        FILTPOS =                    0 / filter wheel position estimated from tm5(raw)
        GAOFSTAT= 'none    '           / status of the CEB gain and offset correction
        BADPXREM= 'off     '           / commanded bad pixel removal on or off
        BADPXDEF=                 4000 / commanded bad pixel default @-8.30 s
        CRREM   = 'off     '           / cosmic ray removal on or off
        CRREMLIM=                 1024 / [1 G. sigma DN] cosmic ray limit @-8.30 s
        GAINHG  =                  256 / global gain corr. high gain @-8.30 s
        GAINLG  =                  256 / global gain corr. how gain  @-8.30 s
        OFFSETHG=                    0 / [DN] global offset corr. high gain @-8.30 s
        OFFSETLG=                    0 / [DN] global offset corr. low gain @-8.30 s
        PRIORITY=                    0 / priority # of image/histogram when downloaded
        SCITABID=                   28 / exposure identifier in sec. science table
        SCITABNR=                    2 / sequential number of SCITABID
        RECSTATE= 'on      '           / recoding on or off
        RECNRBIT=                    8 / bit depth recoding output (sqrt)
        RECLOW  =                    0 / rec. threshold low (clipped 0)
        RECHIGH =                32767 / rec. threshold high(clipped 2^RECNRBIT-1)
        COMBITPP=                    3 / WICOM compression bpp ( COMBITPP*0.04=bpp)
        COMSPLMD= 'provided by user'   / WICOM compression splitb3 mode
        COMSPLVL=                    5 / WICOM compression splitb3 value
        COMWEIMD= 'off     '           / WICOM compression weighting mode
        COMWEIVL= '128,128,128,128,128,128,128,128,128,128' / WICOM sub-band coefficient
        COMSIZE =                 8944 / number of bytes onboard compressed image
        COMSTRIP= 'off     '           / compression type: off=whole image; on=strip
        COMPRESS= 'Lossy-extreme'      / data compression quality (None/Lossless/Lossy)
        COMP_RAT=     1582.71198568873 / compression ratio: uncompressed/compressed size
        EUXCEN  =    378.1416927369173 / [pixel] axis 1 location of solar center in L1
        EUYCEN  =    378.2622041177966 / [pixel] axis 2 location of solar center in L1
        DATAMIN =                    0 / minimum valid physical value
        DATAMAX =                 3486 / maximum valid physical value
        DATAMEAN=    70.92212422688802 / [DN] average pixel value across the image
        RSUN_ARC=    973.9929284938777 / [arcsec] apparent photospheric solar radius
        RSUN_OBS=    973.9929284938777 / [arcsec] apparent photospheric solar radius
        RSUN_REF=            695700000 / [m] assumed physical solar radius
        SOLAR_B0=   -6.677145021348451 / [deg] s/c tilt of solar North pole
        SOLAR_P0=    24.49493954200726 / [deg] s/c celestial North to solar North angle
        SOLAR_EP=    6.020922523528366 / [deg] s/c ecliptic North to solar North angle
        CAR_ROT =    2236.260991931539 / carrington rotation number
        HGLT_OBS=   -6.677145021348451 / [deg] s/c Heliographic latitude (B0 angle)
        HGLN_OBS=    125.2577856554425 / [deg] s/c Heliographic longitude
        CRLT_OBS=   -6.677145021348451 / [deg] s/c Carrington latitude (B0 angle)
        CRLN_OBS=    266.0429046459254 / [deg] s/c Carrington longitude (L0 angle)
        DSUN_OBS=     147330596344.022 / [m] s/c distance from Sun
        DSUN_AU =   0.9848442070373797 / [AU] s/c distance from Sun
        HEEX_OBS=   -85696377695.41205 / [m] s/c Heliocentric Earth Ecliptic X
        HEEY_OBS=    111586055474.3987 / [m] s/c Heliocentric Earth Ecliptic Y
        HEEZ_OBS=     43714845220.3008 / [m] s/c Heliocentric Earth Ecliptic Z
        HCIX_OBS=    30146038635.19564 / [m] s/c Heliocentric Inertial X
        HCIY_OBS=    143192376435.2006 / [m] s/c Heliocentric Inertial Y
        HCIZ_OBS=   -17130799880.04794 / [m] s/c Heliocentric Inertial Z
        HCIX_VOB=   -24672.58385959416 / [m/s] s/c Heliocentric Inertial X Velocity
        HCIY_VOB=    3682.865033100267 / [m/s] s/c Heliocentric Inertial Y Velocity
        HCIZ_VOB=    290.6051210903948 / [m/s] s/c Heliocentric Inertial Z Velocity
        HAEX_OBS=   -132369527358.7171 / [m] s/c Heliocentric Aries Ecliptic X
        HAEY_OBS=    64679546764.96254 / [m] s/c Heliocentric Aries Ecliptic Y
        HAEZ_OBS=    1081238268.853039 / [m] s/c Heliocentric Aries Ecliptic Z
        HEQX_OBS=   -84470625652.15648 / [m] s/c Heliocentric Earth Equatorial X
        HEQY_OBS=    119488717946.4578 / [m] s/c Heliocentric Earth Equatorial Y
        HEQZ_OBS=   -17130799880.04794 / [m] s/c Heliocentric Earth Equatorial Z
        GSEX_OBS=    234594796125.3331 / [m] s/c Geocentric Solar Ecliptic X
        GSEY_OBS=   -111586055474.3987 / [m] s/c Geocentric Solar Ecliptic Y
        GSEZ_OBS=    43714845220.30079 / [m] s/c Geocentric Solar Ecliptic Z
        OBS_VR  =    24947.57547527891 / [m/s] Radial velocity of S/C relative to Sun
        EAR_TDEL=    5.229691555145962 / [s] Time(Sun to Earth) - Time(Sun to S/C)
        SUN_TIME=    491.4419706449787 / [s] Time(Sun to S/C)
        DATE_EAR= '2020-10-21T14:55:18.436' / [UTC] start time of observation Earth
        DATE_SUN= '2020-10-21T14:47:01.764' / [UTC] start time of observation Sun
        INFO_URL= 'http://sidc.be/EUI/data' / Link to additional information
        CHECKSUM= 'Zoh7Zog4Zog4Zog4'   / HDU checksum updated 2020-11-23T16:40:54
        DATASUM = '4177122843'         / data unit checksum updated 2020-11-23T16:40:54
        WAVEUNIT= 'Angstrom'
        BSCALE  =                    1
        BZERO   =                32768
        COMMENT --------- General Description: -----------------------------------------
        COMMENT --------- Instrument and Observation Configuration: --------------------
        COMMENT --------- Description of Data Content: ---------------------------------
        COMMENT --------- Image Relative to Detector and Electronics: ------------------
        COMMENT --------- World Coordinate System Attitude: ----------------------------
        COMMENT --------- Front End Electronics: ---------------------------------------
        COMMENT --------- Temperature: -------------------------------------------------
        COMMENT --------- Telemetry Header: --------------------------------------------
        COMMENT --------- CEB Pixel Preprocessing: -------------------------------------
        COMMENT --------- Data Routing: ------------------------------------------------
        COMMENT --------- Onboard Processing: ------------------------------------------
        COMMENT --------- Derived Image Properties: ------------------------------------
        COMMENT --------- Solar Ephemeris: ---------------------------------------------
        COMMENT --------- Parameters Closing Metadata: ---------------------------------
        HISTORY created by /home/eui/pipeline/telemetry_parser.py --databaseExternal --w
        HISTORY orkingDirectory /tmp/telemetry_parser --configFile /home/eui/config/conf
        HISTORY ig.ini --outputDirectory /data/solo-eui/internal/L0/ --atROBcreated by /
        HISTORY home/eui/pipeline/level0_to_level1.py -i /data/solo-eui/internal/L0/2020
        HISTORY /10/21 --outputDirectory /data/solo-eui/internal/test/L1/ --configFile /
        HISTORY home/eui/config/config.ini --estimatefilter --reprocess
        """)
    header = fits.Header.fromstring(raw_header, sep='\n')
    data = np.random.rand(header['naxis1'], header['naxis2'])
    return Map(data, header)


def test_EUIMap(eui_fsi_map):
    assert isinstance(eui_fsi_map, EUIMap)


def test_is_datasource_for(eui_fsi_map):
    assert eui_fsi_map.is_datasource_for(eui_fsi_map.data, eui_fsi_map.meta)


def test_observer_coordinate(eui_fsi_map):
    obs_coord = eui_fsi_map.observer_coordinate
    assert isinstance(obs_coord, SkyCoord)
    assert obs_coord.obstime.isot == eui_fsi_map.meta['date-avg']


def test_observatory(eui_fsi_map):
    assert eui_fsi_map.observatory == "Solar Orbiter"


def test_measurement(eui_fsi_map):
    assert eui_fsi_map.measurement == u.Quantity(304, 'angstrom')


def test_exposure_time(eui_fsi_map):
    assert eui_fsi_map.exposure_time == u.Quantity(6, 's')


def test_level_number(eui_fsi_map):
    assert eui_fsi_map.processing_level == 1


def test_unit(eui_fsi_map):
    assert eui_fsi_map.unit == u.Unit('ct')


def test_norm_clip(eui_fsi_map):
    # Tests that the default normalizer has clipping disabled
    assert not eui_fsi_map.plot_settings['norm'].clip


def test_wcs(eui_fsi_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    eui_fsi_map.pixel_to_world(0*u.pix, 0*u.pix)
