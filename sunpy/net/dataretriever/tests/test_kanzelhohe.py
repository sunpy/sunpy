#This module was developed with funding provided by
#the Google Summer of Code 2016.
import datetime
import pytest

from astropy import units as u
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument,Level,Wavelength
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.downloader_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
import sunpy.net.dataretriever.sources.kanzelhohe as kanzelhohe

KClient = kanzelhohe.KanzelhoheClient()

@pytest.mark.online
@pytest.mark.parametrize("timerange, wavelength, url_start, url_end",
                         [(TimeRange('2015/01/10 00:00:00', '2015/01/10 12:00:00'), Wavelength(6563*u.AA),
                           'http://cesar.kso.ac.at/halpha2k/recent/2015/kanz_halph_fr_20150110_102629.fts.gz',
                           'http://cesar.kso.ac.at/halpha2k/recent/2015/kanz_halph_fr_20150110_113524.fts.gz')])
def test_get_url_for_timerange(timerange, wavelength, url_start, url_end):
    urls = KClient._get_url_for_timerange(timerange, Wavelength = wavelength )
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[1] == url_end

def test_can_handle_query():
    time = Time('2015/12/30 00:00:00','2015/12/31 00:05:00')
    ans1 = kanzelhohe.KanzelhoheClient._can_handle_query(time, Instrument('kanzelhohe'), Wavelength(6563*u.AA))
    assert ans1 is True
    ans2 = kanzelhohe.KanzelhoheClient._can_handle_query(time, Instrument('swap'))
    assert ans2 is False
    ans3 = kanzelhohe.KanzelhoheClient._can_handle_query(time)
    assert ans3 is False
    ans4 = kanzelhohe.KanzelhoheClient._can_handle_query(time, Instrument('kanzelhohe'), Wavelength(32768*u.AA))
    assert ans4 is True

@pytest.mark.online
def test_query():
    qr = KClient.query(Time('2015/01/10 00:00:00', '2015/01/10 12:00:00'), Wavelength(6563*u.AA))
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 2
    assert qr.time_range()[0] == '2015/01/10'
    assert qr.time_range()[1] == '2015/01/10'

@pytest.mark.online
@pytest.mark.parametrize("time, wavelength",
[(Time('2015/01/02 07:30:00', '2015/01/02 07:38:00'), Wavelength(32768*u.AA))])
def test_get(time, wavelength):
    qr = KClient.query(time, wavelength)
    res = KClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)

@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2016/01/05 07:30:00', '2016/01/05 07:38:00'), a.Instrument('kanzelhohe'), a.Wavelength(5460*u.AA))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
