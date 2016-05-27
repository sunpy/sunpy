#This module was developed with funding provided by
#the Google Summer of Code 2016.
import datetime
import pytest

from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument,Filter
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.downloader_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
import sunpy.net.dataretriever.sources.xrt as xrt

XClient = xrt.XRTClient()

@pytest.mark.online
@pytest.mark.parametrize("timerange,url_start,url_end, filter_",
[(TimeRange('2016/5/18','2016/5/19'),
'http://solar.physics.montana.edu/HINODE/XRT/QL/syn_comp_fits/XRT_Al_mesh_20160518_180137.1.fits',
'http://solar.physics.montana.edu/HINODE/XRT/QL/syn_comp_fits/XRT_Al_mesh_20160518_062007.7.fits' ,'al_mesh')])
def test_get_url_for_timerange(timerange, url_start, url_end, filter_):
    urls = XClient._get_url_for_timerange(timerange, filter = filter_)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end

trange = Time('2015/12/30', '2015/12/31')
def test_can_handle_query():
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('almesh')) is True
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('alpoly')) is True
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('cpoly')) is True
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('tipoly')) is True
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('thin_be')) is True
    assert XClient._can_handle_query(trange, Instrument('xrt'), Filter('mag')) is not True

@pytest.mark.online
def test_query():
    qr = XClient.query(Time('2016/5/18','2016/5/19'), Instrument = 'xrt', filter = 'al_mesh')
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 2
    assert qr.time_range()[0] == '2016/05/18'
    assert qr.time_range()[1] == '2016/05/19'




    
