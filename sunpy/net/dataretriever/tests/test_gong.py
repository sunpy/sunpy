"""
This module tests GONG Client
"""
# This module was developed with funding provided by
# the Google Summer of Code 2016.
import pytest
import datetime

from astropy import units as u
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time.timerange import TimeRange

import sunpy.net.dataretriever.sources.gong as gong

GONGClient = gong.GONGClient()
FClient = gong.FARSIDEClient()

TRANGE = a.Time('2014/6/4 00:00:00', '2014/6/4 00:07:00')


# downloads 4 fits files each of size 3MB.
# Total size = 12MB.
@pytest.mark.online
@pytest.mark.parametrize("time, physobs, instrument, wavelength",
                         [(TRANGE, a.Physobs('INTENSITY'),
                           a.Instrument('maunaloa'), a.Wavelength(656.3*u.nm)),
                          (TRANGE, a.Physobs('INTENSITY'),
                           a.Instrument(''), a.Wavelength(6563*u.AA))])
def test_query(time, physobs, instrument, wavelength):
    qr = GONGClient.search(time, physobs, instrument, wavelength)
    res = GONGClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)


@pytest.mark.parametrize("time, instrument, physobs, wavelength, expected",
                         [(TRANGE, a.Instrument('bigbear'), a.Physobs('INTENSITY'), None, True),
                          (TRANGE, None, a.Physobs('LOS_MAGNETIC_FIELD'), None, True),
                          (TRANGE, a.Instrument('tucson'), None, None, True),
                          (TRANGE, a.Instrument('cerrotololo'), a.Physobs('INTENSITY'),
                           a.Wavelength(6563*u.AA), True),
                          (TRANGE, None, None, None, False)])
def test_can_handle_query(time, instrument, physobs, wavelength, expected):
    assert GONGClient._can_handle_query(time, instrument, physobs, wavelength) is expected


@pytest.mark.online
def test_query_range():
    qr = GONGClient.search(a.Time('2016/6/4 00:00:00', '2016/6/4 00:30:00'),
                           physobs='LOS_MAGNETIC_FIELD', instrument='bigbear')
    assert len(qr) == 3
    assert qr.time_range().start.date() == datetime.date(2016, 6, 4)
    assert qr.time_range().end.date() == datetime.date(2016, 6, 4)


# Downloads 4 fits files each of
# size 1.2MB. Total size = 4.8MB.
@pytest.mark.online
@pytest.mark.parametrize("time, physobs, instrument, wavelength",
                         [(a.Time('2016/6/13 03:00', '2016/6/13 04:00'), a.Physobs('INTENSITY'),
                           a.Instrument('udaipur'), a.Wavelength(676.8*u.nm))])
def test_get(time, physobs, instrument, wavelength):
    qr = GONGClient.search(time, physobs, instrument, wavelength)
    res = GONGClient.get(qr)
    download_list = res.wait()
    assert len(qr) == len(download_list)


# Downloads 2 fits.gz files of size 1.2MB each.
# Total size is 2.5MB.
@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2016/6/4', '2016/6/4 00:10:00'), a.Physobs('INTENSITY'),
                     a.Wavelength(6768*u.AA))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile


@pytest.mark.online
@pytest.mark.parametrize("timerange, url_start, url_end",
                         [(TimeRange('2014/2/1', '2014/2/5'),
                           '201402/mrfqo140201/mrfqo140201t0000.fits',
                           '201402/mrfqo140205/mrfqo140205t0000.fits')])
def test_farside_get_url_for_timerange(timerange, url_start, url_end):
    urls = FClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    domain = 'http://farside.nso.edu/oQR/fqo/'
    assert urls[0] == domain + url_start
    assert urls[-1] == domain + url_end


def test_farside_can_handle_query():
    assert FClient._can_handle_query(a.Time('2015/12/28', '2015/12/30'), a.Instrument('farside'))
    assert FClient._can_handle_query(a.Time('2015/12/28', '2015/12/30'), a.Instrument('Farside'))
    assert not FClient._can_handle_query(a.Time('2015/12/28', '2015/12/30'))
    assert not FClient._can_handle_query(a.Time('2015/12/28', '2015/12/30'), a.Instrument('bbso'))


@pytest.mark.online
def test_farside_query():
    qr = FClient.search(a.Time('2016/1/1', '2016/1/5'), instrument='farside')
    assert isinstance(qr, QueryResponse)
    assert len(qr) == 9
    assert qr.time_range().start.date() == datetime.date(2016, 1, 1)
    assert qr.time_range().end.date() == datetime.date(2016, 1, 5)


# Downloads 7 fits files each of size
# 160KB. Total size ~ 1.2MB
@pytest.mark.online
@pytest.mark.parametrize("time, instrument",
                         [(a.Time('2016/1/1 00:00:00', '2016/1/4 10:00:00'),
                           a.Instrument('farside'))])
def test_farside_get(time, instrument):
    qr = FClient.search(time, instrument)
    res = FClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)


# Downloads 5 fits files each of size 160KB.
# Total size ~ 800KB.
@pytest.mark.online
def test_farside_fido_query():
    qr = Fido.search(a.Time('2016/5/18', '2016/5/20'), a.Instrument('farside'))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
