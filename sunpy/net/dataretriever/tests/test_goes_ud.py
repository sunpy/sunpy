
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
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.sources import goes
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

from hypothesis import given, example
from sunpy.net.tests.strategies import goes_time

GClient = goes.GOESClient()
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
    urls = GClient._get_url_for_timerange(timerange, physobs='IRRADIANCE')
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
        GClient.search(Time("1950/01/01", "1950/02/02"), Instrument('goes'), Physobs('IRRADIANCE'))


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
            GClient.search(time, Instrument('goes'), Physobs('IRRADIANCE'))
    else:
        qr1 = GClient.search(time, Instrument('goes'), Physobs('IRRADIANCE'))
        assert isinstance(qr1, QueryResponse)
        assert qr1.time_range().start == time.start
        assert qr1.time_range().end == time.end


@pytest.mark.online
@pytest.mark.parametrize("time, physobs", [
    (('1983/06/17', '1983/06/18'), 'IRRADIANCE'),
    (('2012/10/04', '2012/10/6'), 'IRRADIANCE'),
])
def test_get(time, physobs):
    qr1 = GClient.search(Time(*time), Instrument('goes'), Physobs(physobs))
    res = GClient.fetch(qr1)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr1)


@pytest.mark.online
def test_new_logic():
    qr = GClient.search(Time('2012/10/4', '2012/10/6'), Instrument('goes'), Physobs('IRRADIANCE'))
    res = GClient.fetch(qr)
    download_list = res.wait(progress=False)
    assert len(download_list) == len(qr)


@pytest.mark.online
@pytest.mark.parametrize("time, physobs",
                         product([('2016/06/15', '2016/06/22'),
                                  ('2013/01/01', '2013/01/08')],
                                   # (UTCnow - datetime.timedelta(days=4), UTCnow)],
                                  ['PARTICLE_FLUX'])
                         )
def test_query_b(time, physobs): # FIXME: Fails
    qr = GClient.search(a.Time(*time), Instrument('goes'), Physobs(physobs))
    assert len(qr) == 8

@pytest.mark.online
@pytest.mark.parametrize("time, physobs, exp",
                         [(Time('2016/06/15', '2016/06/22'), 'IRRADIANCE', 8),
                          (Time(UTCnow, UTCnow), 'IRRADIANCE', 5)
                         ])
def test_query_1(time, physobs, exp):
    qr = GClient.search(time, Instrument('goes'), Physobs(physobs))
    assert len(qr) == exp

#Downloads 2 FITS files with total
#size of ~2.MB
@pytest.mark.online
@pytest.mark.parametrize("time, physobs",
                         [(('2016/06/04', '2016/06/04 00:01:00'), 'INTENSITY')])
def test_query_2(time, physobs):
    qr = GClient.search(Time(*time), Instrument('goes'), Physobs(physobs))
    res = GClient.fetch(qr)
    d_list = res.wait()
    assert len(qr) == len(d_list)


TRANGE = Time(datetime.datetime.now() - datetime.timedelta(days=2), datetime.datetime.now())
@pytest.mark.parametrize("time, physobs, number_of_files",
                         [(TRANGE, Physobs('PARTICLE_FLUX'), 7),
                          (TRANGE, Physobs('IRRADIANCE'), 7),
                          (Time('2016/7/1', '2016/7/3'), Physobs('IRRADIANCE'), 3)])
def test_query_c(time, physobs, number_of_files):
    qr = GClient.search(time, Instrument('goes'), physobs)
    assert len(qr) == number_of_files


TRANGE = Time('2016/6/15', '2016/6/17')
@pytest.mark.parametrize("time, physobs, expected",
                         [(TRANGE, Physobs('PARTICLE_FLUX'), True),
                          (TRANGE, Physobs('INTENSITY'), True),
                          (TRANGE, Physobs('IRRADIANCE'), True),
                          (TRANGE, None, False)])
def test_can_handle_query_b(time, physobs, expected):
    assert GClient._can_handle_query(time, Instrument('goes'), physobs) is expected


@pytest.mark.online
@pytest.mark.parametrize("time, instrument, physobs",
                         [(('2012/10/04', '2012/10/06'), "goes", 'INTENSITY'),
                          (('2013/10/05', '2013/10/07'), "goes", 'INTENSITY')
                         ])
def test_fido(time, instrument, physobs):
    qr = Fido.search(a.Time(*time), a.Instrument(instrument), Physobs(physobs))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
