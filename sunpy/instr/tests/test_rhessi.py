"""
Because rhessi.py is code as well.
"""
import os

import sunpy.map
import sunpy.data.test
import sunpy.instr.rhessi as rhessi

import numpy as np
import pytest

from datetime import datetime

testpath = sunpy.data.test.rootdir


def test_backprojection():
    amap = rhessi.backprojection(os.path.join(testpath, 'hsi_calib_ev_20020220_1106_20020220_1106_25_40.fits'))
    assert isinstance(amap, sunpy.map.GenericMap)
    assert amap.date == datetime(2002, 2, 20, 11, 6, 21)


def test_get_obssumm_dbase_file():
    with pytest.raises(ValueError):
        rhessi.get_obssumm_dbase_file(['2002/01/01', '2002/04/01'])


@pytest.mark.online
def test_get_obssum_filename():
    file_name = rhessi.get_obssum_filename(('2011/04/04', '2011/04/05'))
    # Irregardless of mirror server the osbsum file name should match
    assert file_name[0].split('metadata/catalog/')[1] == 'hsi_obssumm_20110404_042.fits'


@pytest.mark.online
def test_parse_obssum_dbase_file():
    file = rhessi.get_obssumm_dbase_file(('2011/04/04', '2011/04/05'))
    obssum = rhessi.parse_obssumm_dbase_file(file[0])

    assert obssum['filename'][0] == 'hsi_obssumm_20110401_043.fit'
    assert obssum['filename'][-1] == 'hsi_obssumm_20110430_029.fit'

    assert obssum['orb_st'][0] == 0
    assert obssum['orb_st'][-1] == 0

    assert obssum['orb_end'][0] == 0
    assert obssum['orb_end'][-1] == 0

    assert obssum['start_time'][0] == datetime(2011, 4, 1, 0, 0, 0)
    assert obssum['start_time'][-1] == datetime(2011, 4, 30, 0, 0, 0)

    assert obssum['end_time'][0] == datetime(2011, 4, 2, 0, 0, 0)
    assert obssum['end_time'][-1] == datetime(2011, 5, 1, 0, 0, 0)

    assert obssum['status_flag'][0] == 0
    assert obssum['status_flag'][-1] == 0

    assert obssum['npackets'][0] == 0
    assert obssum['npackets'][-1] == 0


@pytest.mark.online
def test_get_parse_obssum_file():
    f = rhessi.get_obssumm_file(('2011/04/04', '2011/04/05'))  # doctest: +SKIP
    header, _data = rhessi.parse_obssumm_file(f[0])
    assert header.get('DATE_OBS') == '2011-04-04T00:00:00.000'
    assert header.get('DATE_END') == '2011-04-05T00:00:00.000'
    assert header.get('TELESCOP') == 'HESSI'


def test_uncompress_countrate():
    # Should only accept bytearr (uncompressed counts must be 0 - 255)
    with pytest.raises(ValueError):
        rhessi.uncompress_countrate(np.array([-1, 300]))

    counts = rhessi.uncompress_countrate(np.array([0, 128, 255]))

    # Valid min, max
    assert counts[0] == 0
    assert counts[2] == 1015792

    # Random test value
    assert counts[1] == 4080
