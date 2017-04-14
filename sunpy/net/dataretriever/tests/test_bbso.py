# This module was developed with funding provided by
# the Google Summer of Code 2016.
import pytest

from sunpy.time import parse_time
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument, Level
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

import sunpy.net.dataretriever.sources.bbso as bbso

BClient = bbso.BBSOClient()


@pytest.mark.online
@pytest.mark.parametrize("timerange,url_start,url_end, level", [(
    TimeRange('2016/4/4 15:28:00', '2016/4/4 16:40:00'),
    'http://www.bbso.njit.edu/pub/archive/2016/04/04/bbso_halph_fl_20160404_152959.fts',
    'http://www.bbso.njit.edu/pub/archive/2016/04/04/bbso_halph_fl_20160404_163016.fts',
    'fl')])
def test_get_url_for_time_range(timerange, url_start, url_end, level):
    urls = BClient._get_url_for_timerange(timerange, level=level)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


trange = Time('2015/12/30', '2015/12/31')


def test_can_handle_query():
    assert BClient._can_handle_query(trange, Instrument('bbso'), Level('fl'))
    assert BClient._can_handle_query(trange, Instrument('bbso'), Level('fr'))
    assert not BClient._can_handle_query(trange,
                                         Instrument('bbso'), Level('mag'))
    assert not BClient._can_handle_query(trange, Level('mag'))


@pytest.mark.online
def test_query():
    qr = BClient.query(
        Time('2016/5/18 15:28:00', '2016/5/18 16:30:00'),
        Instrument='bbso',
        Level='fr')
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 2
    assert qr.time_range().start == parse_time('2016/05/18 15:30:25')
    assert qr.time_range().end == parse_time('2016/05/18 16:00:33')


#This test downloads 2 fits files
#each of size 8.4MB, total size 16.8MB
@pytest.mark.online
@pytest.mark.parametrize("time, instrument, level",
                         [(Time('2016/5/18 15:28:00', '2016/5/18 16:30:00'),
                           Instrument('bbso'), Level('fr'))])
def test_get(time, instrument, level):
    qr = BClient.query(time, instrument, level)
    res = BClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)


#This test downloads 2 fits files
#each of size 8MB, total size 16MB
@pytest.mark.online
def test_fido_query():
    qr = Fido.search(
        a.Time('2016/03/02 17:00:00', '2016/03/02 17:35:00'),
        a.Instrument('bbso'), a.Level('fl'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
