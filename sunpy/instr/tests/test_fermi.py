
import pytest
import collections
from numpy.testing import assert_almost_equal
from sunpy.instr import fermi
from sunpy.time import parse_time



@pytest.mark.online
def test_download_weekly_pointing_file():
    #set a test date
    date = parse_time('2011-10-01')
    file = fermi.download_weekly_pointing_file(date)
    assert file.endswith('FERMI_POINTING_FINAL_174_2011272_2011279_02.fits')


@pytest.mark.online
def test_detector_angles():
    #set a test date
    date = parse_time('2012-02-15')
    file = fermi.download_weekly_pointing_file(date)
    det=fermi.get_detector_sun_angles_for_date(date,file)
    assert len(det) == 13
    #assert type(det) == collections.OrderedDict
    assert_almost_equal(det['n0'][0].value, 20.30309,decimal=1)
    assert_almost_equal(det['n1'][0].value, 30.30430, decimal=1)
    assert_almost_equal(det['n2'][0].value, 74.86032, decimal=1)
    assert_almost_equal(det['n3'][0].value, 31.24400, decimal=1)
    assert_almost_equal(det['n4'][0].value, 75.10403, decimal=1)
    assert_almost_equal(det['n5'][0].value, 60.40967, decimal=1)
    assert_almost_equal(det['n6'][0].value, 46.14087, decimal=1)
    assert_almost_equal(det['n7'][0].value, 69.71780, decimal=1)
    assert_almost_equal(det['n8'][0].value, 106.08064, decimal=1)
    assert_almost_equal(det['n9'][0].value, 68.543067, decimal=1)
    assert_almost_equal(det['n10'][0].value, 105.76825, decimal=1)
    assert_almost_equal(det['n11'][0].value, 119.69057, decimal=1)

    det2=fermi.get_detector_sun_angles_for_time(parse_time('2012-02-15 02:00'),file)
    assert len(det2) == 13
    assert type(det2) == dict
    assert_almost_equal(det2['n0'].value, 87.24744,decimal=1)
    assert_almost_equal(det2['n1'].value, 69.90883,decimal=1)
    assert_almost_equal(det2['n10'].value, 123.56429,decimal=1)
    assert_almost_equal(det2['n11'].value, 167.26615,decimal=1)
    assert_almost_equal(det2['n2'].value, 59.82642,decimal=1)
    assert_almost_equal(det2['n3'].value, 69.18959,decimal=1)
    assert_almost_equal(det2['n4'].value, 56.83158,decimal=1)
    assert_almost_equal(det2['n5'].value, 12.49959,decimal=1)
    assert_almost_equal(det2['n6'].value, 115.31259,decimal=1)
    assert_almost_equal(det2['n7'].value, 129.49283,decimal=1)
    assert_almost_equal(det2['n8'].value, 121.91083,decimal=1)
    assert_almost_equal(det2['n9'].value, 130.04144,decimal=1)
        
