"""
This module tests XRT Client.
"""
#This module was developed with funding provided by
#the Google Summer of Code 2016.
import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument, Filter
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.downloader_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
import sunpy.net.dataretriever.sources.xrt as xrt

XClient = xrt.XRTClient()

@pytest.mark.online
@pytest.mark.parametrize("timerange,url_start,url_end, filter_",
[(TimeRange('2016/5/18', '2016/5/19'),
'http://solar.physics.montana.edu/HINODE/XRT/QL/syn_comp_fits/XRT_Al_mesh_20160518_180137.1.fits',
'http://solar.physics.montana.edu/HINODE/XRT/QL/syn_comp_fits/XRT_Al_mesh_20160518_062007.7.fits', 'al_mesh')])
def test_get_url_for_timerange(timerange, url_start, url_end, filter_):
    urls = XClient._get_url_for_timerange(timerange, filter=filter_)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

trange = Time('2015/12/30', '2015/12/31')
def test_can_handle_query():
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('almesh'))
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('alpoly'))
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('cpoly'))
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('tipoly'))
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('thin_be'))
    assert not XClient._can_handle_query(trange, Instrument('xrt'), Filter('mag'))

@pytest.mark.online
def test_query():
    qr = XClient.query(Time('2016/5/18', '2016/5/19'), Instrument='xrt', filter='al_mesh')
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 2
    assert qr.time_range()[0] == '2016/05/18'
    assert qr.time_range()[1] == '2016/05/19'

@pytest.mark.online
@pytest.mark.parametrize("time, instrument, filter_",
[(Time('2016/5/18', '2016/5/19'), Instrument('xrt'), Filter('al_mesh'))])
def test_get(time, instrument, filter_):
    qr = XClient.query(time, instrument, filter_)
    res = XClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)

#Downloads two fits files each of size 4MB
#Total size = 8MB
@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2016/5/18', '2016/5/19'), a.Instrument('xrt'), a.Filter('al_mesh'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
