
"""
This module tests GOES Client.
"""
#This module was developed by funding
#provided by Google Summer of Code 2016.
import pytest
import datetime
from itertools import product

from sunpy.time.timerange import TimeRange, parse_time
from sunpy.net.vso.attrs import Time, Instrument, Physobs
import sunpy.net.dataretriever.sources.goes as goes
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

from hypothesis import given, example
from sunpy.net.tests.strategies import goes_time

LCClient = goes.GOESClient()
UTCnow = datetime.datetime.utcnow()

@pytest.mark.parametrize(
    "timerange,url_start,url_end",
    [(TimeRange('1995/06/03', '1995/06/05'),
      'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950603.fits',
      'http://umbra.nascom.nasa.gov/goes/fits/1995/go07950605.fits'),
     (TimeRange('2008/06/02', '2008/06/04'),
      'http://umbra.nascom.nasa.gov/goes/fits/2008/go1020080602.fits',
      'http://umbra.nascom.nasa.gov/goes/fits/2008/go1020080604.fits')])
def test_get_url_for_time_range(timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange, physobs='IRRADIANCE')
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


@given(goes_time())
def test_can_handle_query(time):
    ans1 = goes.GOESClient._can_handle_query(time, Instrument('goes'), Physobs('IRRADIANCE'))
    assert ans1 is True
    ans2 = goes.GOESClient._can_handle_query(time)
    assert ans2 is False
    ans3 = goes.GOESClient._can_handle_query(time, Instrument('eve'))
    assert ans3 is False

def test_no_satellite():
    with pytest.raises(ValueError):
        LCClient.search(Time("1950/01/01", "1950/02/02"), Instrument('goes'), Physobs('IRRADIANCE'))

@example(a.Time("2006-08-01", "2006-08-01"))
@example(a.Time("1980-01-03", "1980-01-05"))
# This example tests a time range with a satellite jump and no overlap
@example(a.Time("2009-11-30", "2009-12-3"))
@given(goes_time())
@pytest.mark.online
def test_query_a(time):
    tr = TimeRange(time.start, time.end)
    if parse_time("1980-01-03") in tr:
        with pytest.raises(ValueError):
            LCClient.search(time, Instrument('goes'), Physobs('IRRADIANCE'))
    else:
        qr1 = LCClient.search(time, Instrument('goes'), Physobs('IRRADIANCE'))
        assert isinstance(qr1, QueryResponse)
        assert qr1.time_range().start == time.start
        assert qr1.time_range().end == time.end

@pytest.mark.online
@pytest.mark.parametrize("time, instrument, physobs", [
    (('1983/06/17', '1983/06/18'), 'goes', 'IRRADIANCE'),
    (('2012/10/04', '2012/10/6'), 'goes', 'IRRADIANCE'),
])
def test_get(time, instrument, physobs):
    qr1 = LCClient.search(Time(*time), Instrument(instrument), Physobs(physobs))
    res = LCClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@pytest.mark.online

def test_new_logic():
    qr = LCClient.search(Time('2012/10/4', '2012/10/6'), Instrument('goes'), Physobs('IRRADIANCE'))
    res = LCClient.fetch(qr)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr)
@pytest.mark.online
@pytest.mark.parametrize("time, instrument, physobs",
                         product([('2016/06/15', '2016/06/22'),
                                   (UTCnow - datetime.timedelta(days=7), UTCnow)],
                                  'goes', 'PARTICLE_FLUX')
                         )
def test_query_b(time, instrument, physobs): # FIXME: Fails
    qr = LCClient.query(a.Time(*time), Instrument(instrument), Physobs(physobs))
    assert len(qr) == 8

@pytest.mark.online
@pytest.mark.parametrize("time, instrument, physobs, exp",
                         [(Time('2016/06/15', '2016/06/22'), 'goes', 'IRRADIANCE', 8),
                          (Time(UTCnow, UTCnow), 'goes', 'IRRADIANCE', 4) # FIXME: it gives 5
                         ])
def test_query_1(time, instrument, physobs, exp):
    qr = LCClient.query(time, Instrument(instrument), Physobs(physobs))
    assert len(qr) == exp

#Downloads 2 FITS files with total
#size of ~2.MB
@pytest.mark.online
@pytest.mark.parametrize("time, instrument, physobs",
                         [(('2016/06/04', '2016/06/04 00:01:00'), 'goes', 'INTENSITY')])
def test_query_2(time, instrument, physobs):
    qr = LCClient.query(Time(*time), Instrument(instrument), Physobs(physobs))
    res = LCClient.get(qr)
    d_list = res.wait()
    assert len(qr) == len(d_list)

@pytest.mark.parametrize("time, instrument, physobs, number_of_files",
[(Time(datetime.datetime.now() - datetime.timedelta(days=2), datetime.datetime.now()), Instrument('goes'),
  Physobs('PARTICLE_FLUX'), 8),
 (Time(datetime.datetime.now() - datetime.timedelta(days=2), datetime.datetime.now()), Instrument('goes'),
  Physobs('IRRADIANCE'), 16),
 (Time('2016/7/1', '2016/7/3'), Instrument('goes'), Physobs('IRRADIANCE'), 10)])
def test_query_c(time, instrument, physobs, number_of_files):
    qr = GClient.query(time, instrument, physobs)
    assert len(qr) == number_of_files

TRANGE = Time('2016/6/15', '2016/6/17')
@pytest.mark.parametrize("time, instrument, physobs, expected",
                         [(TRANGE, Instrument('goes'), Physobs('PARTICLE_FLUX'), True),
                          (TRANGE, Instrument('goes'), Physobs('INTENSITY'), True),
                          (TRANGE, Instrument('goes'), Physobs('IRRADIANCE'), True),
                          (TRANGE, Instrument('goes'), None, False)])
def test_can_handle_query_b(time, instrument, physobs, expected):
    assert GClient._can_handle_query(time, instrument, physobs) is expected

@pytest.mark.online
@pytest.mark.parametrize("time, instrument, physobs",
                         [(("2012/10/04", "2012/10/06"), "goes", 'INTENSITY'),
                          (('2013/10/05', '2013/10/07'), "goes", 'INTENSITY')
                         ])
def test_fido(time, instrument, physobs):
    print(a.Time(*time))
    qr = Fido.search(a.Time(*time), a.Instrument(instrument), Physobs(physobs))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
