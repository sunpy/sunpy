"""
This module tests the Kanzelhohe Client
"""
#This module was developed with funding provided by
#the Google Summer of Code 2016.
import pytest

from astropy import units as u
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time, Instrument, Wavelength
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
    urls = KClient._get_url_for_timerange(timerange, wavelength = wavelength)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[1] == url_end
TRANGE = Time('2015/12/30 00:00:00', '2015/12/31 00:05:00')
@pytest.mark.parametrize("time, instrument, wavelength, expected",
                         [(TRANGE, Instrument('kanzelhohe'), Wavelength(6563*u.AA), True),
                          (TRANGE, Instrument('swap'), None, False),
                          (TRANGE, None, None, False),
                          (TRANGE, Instrument('kanzelhohe'), Wavelength(32768*u.AA), True)])
def test_can_handle_query(time, instrument, wavelength, expected):
    assert KClient._can_handle_query(time, instrument, wavelength) is expected

@pytest.mark.online
def test_query():
    qr = KClient.query(Time('2015/01/10 00:00:00', '2015/01/10 12:00:00'), Wavelength(6563*u.AA))
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 2
    assert qr.time_range()[0] == '2015/01/10'
    assert qr.time_range()[1] == '2015/01/10'

@pytest.mark.online
@pytest.mark.parametrize("time, wavelength",
[(Time('2015/01/02 07:30:00', '2015/01/02 07:38:00'), Wavelength(3276.8*u.nm))])
#This test downloads 3 files
#Each file is 4.5MB, total size
#is 13.4MB
def test_get(time, wavelength):
    qr = KClient.query(time, wavelength)
    res = KClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)

This test downloads 3 files
Each file is 4.5MB, total size
is 13.4MB
@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2016/01/05 07:30:00', '2016/01/05 07:38:00'), a.Instrument('kanzelhohe'), a.Wavelength(546.0*u.nm))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
