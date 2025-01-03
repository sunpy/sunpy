"""Tests for PSP/WISPR"""
from textwrap import dedent

import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

from sunpy.map import Map
from sunpy.map.sources import WISPRMap
from .helpers import _test_private_date_setters


@pytest.fixture
def wispr_map(scope="module"):
    # Please do not edit this header, it is a 1:1 copy of a SOLO EUI FSI 304 map
    raw_header = dedent("""\
        SIMPLE  =                    T / Written by IDL:  Wed Sep  9 14:48:00 2020
        BITPIX  =                   32 /  32-bit twos complement binary integer
        NAXIS   =                    2 /
        NAXIS1  =                  960 /
        NAXIS2  =                 1024 /
        FILENAME= 'psp_L1_wispr_20200125T000229_V1_2302.fits' /
        FILE_RAW= 'CmpB.12EE4DE3'      /
        APID    = '41b     '           /
        DATE    = '2020-09-09T18:47:59.761' /
        DATE-OBS= '2020-01-25T00:02:29.618' /
        DATE-BEG= '2020-01-25T00:02:29.618' /
        DATE-AVG= '2020-01-25T00:08:20.842' /
        DATE-END= '2020-01-25T00:14:12.067' /
        OBT_BEG =        317606692.959 /
        OBT_END =        317607252.950 /
        NSUMEXP =                    5 /
        NSUMBAD =                    0 /
        XPOSURE =              700.000 /
        IP_TIMET=                  139 /
        READTIME=              2.45760 /
        TELAPSE =              702.449 /
        TIMESYS = 'UTC     '           /
        LEVEL   = 'L1      '           /
        CREATOR = 'wispr_reduce.pro,v 1.63' /
        ORIGIN  = 'NRL     '           /
        DETECTOR=                    2 /
        REGION  =                    2 /
        READOUT0=                    3 /
        CAMERA  = 'FMCFMD  '           /
        CICUCODE= '129     '           /
        DRBUCODE= '130     '           /
        DRB1UCOD= '130     '           /
        DRB2UCOD= '130     '           /
        OBJECT  = 'OuterFFV'           /
        OBS_MODE= 'SYNOPTIC'           /
        TARGET  = 'VexpOuterL.os'      /
        TIMELINE= '584,0x1a050'        /
        STUDY_ID=                  742 /
        STUDYCNT=                    1 /
        VERS_CAL= '5; 57   '           /
        PRIORITY=                    1 /
        VERSION =                    1 /
        FINAL   =                    T /
        BSCALE  =              1.00000 /
        BZERO   =              0.00000 /
        BUNIT   = 'DN      '           /
        DATAMIN =              856.000 /  Minimum Value Not Equal to Zero before BSCALE
        DATAMAX =              249069. /  Maximum Value before BSCALE
        BLANK   =                61166 /
        PXBEG1  =                    1 /
        PXEND1  =                 2048 /
        PXBEG2  =                    1 /
        PXEND2  =                 1920 /
        R1COL   =                    1 /
        R2COL   =                 1920 /
        R1ROW   =                    1 /
        R2ROW   =                 2048 /
        RECTIFY =                    T /
        RECTROTA=                    6 /
        NBIN1   =                    2 /
        NBIN2   =                    2 /
        NBIN    =              4.00000 /
        COMPRESS= 'Lossless'           /
        COMP_RAT=             0.675516 /
        COSMICR =                    T /
        DSTART1 =                    1 /
        DSTOP1  =                  955 /
        DSTART2 =                    1 /
        DSTOP2  =                 1019 /
        IMGCTR  =                 8151 /
        CCIMGSEQ=                 2505 /
        CCEXPCTL=                  999 /
        CCIMGCTR=                   15 /
        COSMICS =                35210 /
        CAMSATPX=                    0 /
        CAMSTRVX=                18638 /
        IP_INIT =                    2 /
        IPIMBUFF=                30840 /
        IPCMBUFF=                  320 /
        IPSATPX =                16000 /
        IPSTRVPX=                  300 /
        IPMINPX =                    0 /
        IPMAXPX =                16383 /
        IPTRUNC =                    0 /
        IPBIAS  =                    2 /
        IPCRSMUL=              4.00000 /
        IPBITSPP=                   19 /
        IPMASK  =                    0 /
        IPMASKCR= '144,879,0,959'      /
        LEDSTATE= 'Off     '           /
        LEDDAC  =                    0 /
        OFFSET  =                 2256 /
        GAINCMD =                   12 /
        GAINMODE= 'HIGH    '           /
        IP_FUNC =                16584 /
        IPCMPCTL=                  305 /
        XFBYTES =              1577140 /
        SCFLAGS = '0x0015ff3e'         /
        CIE_T   =             -39.9500 /
        DRB1_T  =             -39.3083 /
        DET1_T  =             -70.9754 /
        DRB2_T  =             -39.0940 /
        DET2_T  =             -70.9950 /
        ISPREG0 =                   28 /
        ISPREG1 =                61166 /
        ISPREG2 =                61166 /
        ISPREG3 =                61166 /
        ISPREG4 =                61166 /
        ISPREG5 =                61166 /
        ISPREG6 =                  825 /
        ISPREG7 =                  742 /
        UCODREG0=                61166 /
        UCODREG1=                  128 /
        UCODREG2=                    0 /
        UCODREG3=                    0 /
        UCODREG4=                 1920 /
        UCODREG5=                61166 /
        UCODREG6=                61166 /
        UCODREG7=                61166 /
        DSUN_OBS=        31678123000.0 /
        SC_YAW  =            0.0127380 /
        SC_PITCH=         -0.000853018 /
        SC_ROLL =             -3.83791 /
        DATAZER =                    0 /  Number of Zero Pixels
        DATASAT =                    0 /  Number of Saturated Pixels
        DSATVAL =               524287 /  Value used as saturated
        DATAAVG =              107440. /  Between 0 and DSATVAL
        DATAMDN =              96257.0 /  Between 0 and DSATVAL
        DATASIG =              43882.8 /  Between 0 and DSATVAL
        DATAP01 =                 5409 /
        DATAP10 =                63754 /
        DATAP25 =                77544 /
        DATAP50 =                96257 /
        DATAP75 =               128137 /
        DATAP90 =               171321 /
        DATAP95 =               200730 /
        DATAP98 =               230366 /
        DATAP99 =               237914 /
        OBSRVTRY= 'Parker Solar Probe' /
        INSTRUME= 'WISPR   '           /
        WCSNAME = 'Helioprojective Zenith Polynomial' /
        CTYPE1  = 'HPLN-ZPN'           /
        CTYPE2  = 'HPLT-ZPN'           /
        CUNIT1  = 'deg     '           /
        CUNIT2  = 'deg     '           /
        CRPIX1  =              492.333 /
        CRPIX2  =              513.174 /
        PC1_1   =       0.998328091852 /
        PC1_2   =      0.0578015658869 /
        PC2_1   =     -0.0578015658869 /
        PC2_2   =       0.998328091852 /
        CDELT1  =            0.0570821 /
        CDELT2  =            0.0570821 /
        CRVAL1  =        76.6509814623 /
        CRVAL2  =       -13.2570275339 /
        ATT_FILE= 'spp_2020_025_02.ah.bc' /
        PV1_1   =        0.00000000000 /
        PV1_2   =        90.0000000000 /
        PV1_3   =        180.000000000 /
        PV2_0   =   -0.000429731007898 /
        PV2_1   =        1.01042997837 /
        PV2_2   =     -0.0547256991267 /
        PV2_3   =      0.0849808976054 /
        PV2_4   =      -0.175265997648 /
        PV2_5   =      0.0547404997051 /
        LONPOLE =        180.000000000 /
        LATPOLE =        0.00000000000 /
        CTYPE1A = 'RA---ZPN'           /
        CTYPE2A = 'DEC--ZPN'           /
        CUNIT1A = 'deg     '           /
        CUNIT2A = 'deg     '           /
        PC1_1A  =       0.932533562542 /
        PC1_2A  =      -0.361083307193 /
        PC2_1A  =       0.361083307193 /
        PC2_2A  =       0.932533562542 /
        CRPIX1A =              492.333 /
        CRPIX2A =              513.174 /
        CRVAL1A =        170.519696093 /
        CRVAL2A =       -2.59678257890 /
        CDELT1A =           -0.0570821 /
        CDELT2A =            0.0570821 /
        PV1_1A  =        0.00000000000 /
        PV1_2A  =        90.0000000000 /
        PV1_3A  =        180.000000000 /
        PV2_0A  =   -0.000429731007898 /
        PV2_1A  =        1.01042997837 /
        PV2_2A  =     -0.0547256991267 /
        PV2_3A  =      0.0849808976054 /
        PV2_4A  =      -0.175265997648 /
        PV2_5A  =      0.0547404997051 /
        LONPOLEA=        180.000000000 /
        LATPOLEA=        0.00000000000 /
        RSUN_ARC=        4524.67854698 /
        RSUN_REF=        695507968.000 /
        SOLAR_EP=        15.9285776205 /
        CAR_ROT =                 2226 /
        HGLT_OBS=       0.390715276104 /
        HGLN_OBS=       -54.7633009507 /
        CRLT_OBS=       0.390715276104 /
        CRLN_OBS=        59.8905844055 /
        HEEX_OBS=        18187298676.0 /
        HEEY_OBS=       -25963488915.7 /
        HEEZ_OBS=       -241652330.942 /
        HCIX_OBS=        31492170156.9 /
        HCIY_OBS=       -3624289322.90 /
        HCIZ_OBS=        216174389.143 /
        HCIX_VOB=       -38361.9314129 /
        HCIY_VOB=        71493.5202110 /
        HCIZ_VOB=       -4786.08742968 /
        HAEX_OBS=        11254842547.2 /
        HAEY_OBS=        29634581727.6 /
        HAEZ_OBS=       -243044691.612 /
        HEQX_OBS=        18289512873.5 /
        HEQY_OBS=       -25891812859.2 /
        HEQZ_OBS=        216174389.143 /
        EAR_TIME=        383.108719220 /
        SUN_TIME=        105.742395411 /
        HISTORY make_wispr_hdr.pro,v 1.57
        HISTORY wincs-gse.nrl.navy.mil
        HISTORY FSW version 932
        HISTORY wispr_reduce.pro,v 1.63
        HISTORY wispr_readrawsci.pro,v 1.17
        HISTORY Id: get_wispr_pointing.pro,v 1.29 2020/06/02 19:29:24 nathan Exp
        """)
    header = fits.Header.fromstring(raw_header, sep='\n')
    data = np.random.rand(header['naxis1'], header['naxis2'])
    return Map(data, header)


def test_WISPRMap(wispr_map):
    assert isinstance(wispr_map, WISPRMap)


def test_is_datasource_for(wispr_map):
    assert wispr_map.is_datasource_for(wispr_map.data, wispr_map.meta)


def test_observer_coordinate(wispr_map):
    obs_coord = wispr_map.observer_coordinate
    assert isinstance(obs_coord, SkyCoord)
    assert obs_coord.obstime.isot == wispr_map.meta['date-avg']


def test_observatory(wispr_map):
    assert wispr_map.observatory == "Parker Solar Probe"


def test_reference_date(wispr_map):
    assert wispr_map.reference_date.isot == "2020-01-25T00:08:20.842"


def test_date(wispr_map):
    assert wispr_map.date.isot == "2020-01-25T00:02:29.618"


def test_private_date_setters(wispr_map):
    _test_private_date_setters(wispr_map)


def test_measurement(wispr_map):
    assert wispr_map.measurement is None


def test_wavelength(wispr_map):
    assert wispr_map.wavelength is None


def test_exposure_time(wispr_map):
    assert wispr_map.exposure_time == u.Quantity(700, 's')


def test_processing_level(wispr_map):
    assert wispr_map.processing_level == 1

    for value, expected in [
            ('L1', 1),
            ('L2', 2),
            ('L2b', '2b'),
            ('L3', 3),
            ('LW', 'W')]:
        wispr_map.meta['level'] = value
        assert wispr_map.processing_level == expected


def test_detector(wispr_map):
    assert wispr_map.detector == 'Outer'

    wispr_map.meta['DETECTOR'] = 1
    assert wispr_map.detector == 'Inner'

    wispr_map.meta['DETECTOR'] = 'other_val'
    assert wispr_map.detector == 'other_val'


def test_unit(wispr_map):
    assert wispr_map.unit == u.Unit('DN')


def test_norm_clip(wispr_map):
    # Tests that the default normalizer has clipping disabled
    assert not wispr_map.plot_settings['norm'].clip


def test_name(wispr_map):
    assert wispr_map.name == 'WISPR Outer 2020-01-25 00:02:29'


def test_wcs(wispr_map):
    # Smoke test that WCS is valid and can transform from pixels to world coordinates
    wispr_map.pixel_to_world(0*u.pix, 0*u.pix)
