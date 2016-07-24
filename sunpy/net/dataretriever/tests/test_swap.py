"""
This module tests SWAPClient
"""
#This module was developed with funding provided by
#the Google Summer of Code 2016.
import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument, Level
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.downloader_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

import sunpy.net.dataretriever.sources.swap as swap

LCClient = swap.SWAPClient()

@pytest.mark.online
@pytest.mark.parametrize("timerange,url_start,url_end, level",
                         [(TimeRange('2015/12/30 00:00:00','2015/12/30 23:59:59'),
                           'http://proba2.oma.be/swap/data/bsd/2015/12/30/swap_lv1_20151230_000044.fits',
                           'http://proba2.oma.be/swap/data/bsd/2015/12/30/swap_lv1_20151230_235935.fits', 1)])
def test_get_url_for_time_range(timerange, url_start, url_end, level):
    urls = LCClient._get_url_for_timerange(timerange, level = level)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

TRANGE = Time('2015/12/30 00:00:00', '2015/12/31 00:05:00')
@pytest.mark.parametrize("time, instrument, level, expected",
                         [(TRANGE, Instrument('swap'), Level(1), True),
                          (TRANGE, Instrument('swap'), Level('q'), True),
                          (TRANGE, Instrument('swap'), Level('s'), False),
                          (TRANGE, Instrument('swap'), None, False),
                          (TRANGE, None, None, False),
                          (TRANGE, Instrument('eve'), None, False)])
def test_can_handle_query(time, instrument, level, expected):
    assert LCClient._can_handle_query(time, instrument, level) is expected
##    trange = Time('2015/12/30 00:00:00', '2015/12/31 00:05:00')
##    assert swap.SWAPClient._can_handle_query(trange, Instrument('swap'), a.Level(1))
##    assert swap.SWAPClient._can_handle_query(trange, Instrument('swap'), a.Level('q'))
##    assert not swap.SWAPClient._can_handle_query(trange, Instrument('swap'), a.Level('s'))
##    assert not swap.SWAPClient._can_handle_query(trange, Instrument('swap'))
##    assert not swap.SWAPClient._can_handle_query(trange)
##    assert not swap.SWAPClient._can_handle_query(trange, Instrument('eve'))

@pytest.mark.online
def test_query():
    qr1 = LCClient.query(Time('2015-12-30 00:00:00', '2015-12-30 00:05:00'), Level = 1)
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 4
    assert qr1.time_range()[0] == '2015/12/30'
    assert qr1.time_range()[1] == '2015/12/30'

@pytest.mark.online
@pytest.mark.parametrize("time, instrument, level",
[(Time('2015/12/30 00:00:00', '2015/12/30 00:05:00'), Instrument('swap'), Level('q'))])
#This test will download 4 JP2 files
#each of size 170KB.
def test_get(time, instrument, level):
    qr1 = LCClient.query(time, instrument, level)
    res = LCClient.get(qr1)
    download_list = res.wait()
    assert len(download_list) == len(qr1)

#This test will download 2 fits files
#each of size 2MB
@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2015/12/28 00:00:00', '2015/12/28 00:03:00'), a.Instrument('swap'), a.Level(1))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
