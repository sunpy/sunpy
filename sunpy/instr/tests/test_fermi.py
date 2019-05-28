import pytest
from numpy.testing import assert_almost_equal

import astropy.units as u

from sunpy.instr import fermi
from sunpy.time import parse_time


@pytest.mark.remote_data
def test_download_weekly_pointing_file():
    # set a test date
    date = parse_time('2011-10-01')
    afile = fermi.download_weekly_pointing_file(date)
    assert isinstance(afile, str)
    assert afile.endswith('.fits')


@pytest.mark.remote_data
@pytest.mark.flaky(reruns=5)
def test_detector_angles():
    # set a test date
    date = parse_time('2012-02-15')
    file = fermi.download_weekly_pointing_file(date)
    det = fermi.get_detector_sun_angles_for_date(date, file)
    assert len(det) == 13
    assert_almost_equal(det['n0'][0].value, 21.79, decimal=1)
    assert_almost_equal(det['n1'][0].value, 30.45, decimal=1)
    assert_almost_equal(det['n2'][0].value, 74.44, decimal=1)
    assert_almost_equal(det['n3'][0].value, 30.58, decimal=1)
    assert_almost_equal(det['n4'][0].value, 73.93, decimal=1)
    assert_almost_equal(det['n5'][0].value, 58.78, decimal=1)
    assert_almost_equal(det['n6'][0].value, 47.56, decimal=1)
    assert_almost_equal(det['n7'][0].value, 70.89, decimal=1)
    assert_almost_equal(det['n8'][0].value, 106.54, decimal=1)
    assert_almost_equal(det['n9'][0].value, 70.17, decimal=1)
    assert_almost_equal(det['n10'][0].value, 106.95, decimal=1)
    assert_almost_equal(det['n11'][0].value, 121.32, decimal=1)

    det2 = fermi.get_detector_sun_angles_for_time(
        parse_time('2012-02-15 02:00'), file)
    assert len(det2) == 13
    assert type(det2) == dict
    assert_almost_equal(det2['n0'].value, 83.54, decimal=1)
    assert_almost_equal(det2['n1'].value, 66.50, decimal=1)
    assert_almost_equal(det2['n10'].value, 123.39, decimal=1)
    assert_almost_equal(det2['n11'].value, 170.89, decimal=1)
    assert_almost_equal(det2['n2'].value, 58.84, decimal=1)
    assert_almost_equal(det2['n3'].value, 66.44, decimal=1)
    assert_almost_equal(det2['n4'].value, 57.06, decimal=1)
    assert_almost_equal(det2['n5'].value, 8.85, decimal=1)
    assert_almost_equal(det2['n6'].value, 111.95, decimal=1)
    assert_almost_equal(det2['n7'].value, 127.12, decimal=1)
    assert_almost_equal(det2['n8'].value, 122.93, decimal=1)
    assert_almost_equal(det2['n9'].value, 126.82, decimal=1)


def test_met_to_utc():
    time = fermi.met_to_utc(500000000)
    assert (time - parse_time('2016-11-05T00:53:16.000')) < 1e-7 * u.s
