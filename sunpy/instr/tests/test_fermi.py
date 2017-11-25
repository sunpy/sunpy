import pytest
from numpy.testing import assert_almost_equal
from sunpy.instr import fermi
from sunpy.time import parse_time
from sunpy.extern import six


@pytest.mark.remote_data
def test_download_weekly_pointing_file():
    # set a test date
    date = parse_time('2011-10-01')
    afile = fermi.download_weekly_pointing_file(date)
    assert isinstance(afile, six.string_types)
    assert afile.endswith('.fits')


@pytest.mark.remote_data
def test_detector_angles():
    # set a test date
    date = parse_time('2012-02-15')
    file = fermi.download_weekly_pointing_file(date)
    det = fermi.get_detector_sun_angles_for_date(date, file)
    assert len(det) == 13
    assert_almost_equal(det['n0'][0].value, 21.73944, decimal=1)
    assert_almost_equal(det['n1'][0].value, 30.62983, decimal=1)
    assert_almost_equal(det['n2'][0].value, 74.67486, decimal=1)
    assert_almost_equal(det['n3'][0].value, 30.46062, decimal=1)
    assert_almost_equal(det['n4'][0].value, 73.89734, decimal=1)
    assert_almost_equal(det['n5'][0].value, 58.99893, decimal=1)
    assert_almost_equal(det['n6'][0].value, 47.31091, decimal=1)
    assert_almost_equal(det['n7'][0].value, 70.63391, decimal=1)
    assert_almost_equal(det['n8'][0].value, 106.30992, decimal=1)
    assert_almost_equal(det['n9'][0].value, 70.07033, decimal=1)
    assert_almost_equal(det['n10'][0].value, 106.97884, decimal=1)
    assert_almost_equal(det['n11'][0].value, 121.09603, decimal=1)

    det2 = fermi.get_detector_sun_angles_for_time(
        parse_time('2012-02-15 02:00'), file)
    assert len(det2) == 13
    assert type(det2) == dict
    assert_almost_equal(det2['n0'].value, 83.76092, decimal=1)
    assert_almost_equal(det2['n1'].value, 66.65847, decimal=1)
    assert_almost_equal(det2['n10'].value, 123.28952, decimal=1)
    assert_almost_equal(det2['n11'].value, 170.69869, decimal=1)
    assert_almost_equal(det2['n2'].value, 58.78532, decimal=1)
    assert_almost_equal(det2['n3'].value, 66.69068, decimal=1)
    assert_almost_equal(det2['n4'].value, 57.16402, decimal=1)
    assert_almost_equal(det2['n5'].value, 9.04924, decimal=1)
    assert_almost_equal(det2['n6'].value, 112.21230, decimal=1)
    assert_almost_equal(det2['n7'].value, 127.35783, decimal=1)
    assert_almost_equal(det2['n8'].value, 122.98894, decimal=1)
    assert_almost_equal(det2['n9'].value, 126.95987, decimal=1)
