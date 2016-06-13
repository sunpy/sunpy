#This module was developed with funding provided by
#the Google Summer of Code 2016.
import datetime
import pytest

from astropy import units as u
from sunpy.time.timerange import TimeRange
from sunpy.net.vso.attrs import Time,Instrument,Physobs, Wavelength
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.downloader_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a

import sunpy.net.dataretriever.sources.gong as gong

GONGClient = gong.GONGClient()

##downloads 4 fits files each of size 3MB.
##Total size = 12MB.
trange = Time('2014/6/4 00:00:00', '2014/6/4 00:07:00')
@pytest.mark.online
@pytest.mark.parametrize("time, physobs, instrument, wavelength",
                         [(trange, Physobs('INTENSITY'), Instrument('M'), Wavelength(6563*u.AA)),
                          (trange, Physobs('INTENSITY'), Instrument(''), Wavelength(6563*u.AA))])
def test_query(time, physobs, instrument, wavelength):
    qr = GONGClient.query(time, physobs, instrument, wavelength)
    res = GONGClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)

def test_can_handle_query():
    assert GONGClient._can_handle_query(trange, Instrument('bb'), Physobs('INTENSITY'))
    assert GONGClient._can_handle_query(trange, Physobs('LOS_MAGNETIC_FIELD'))
    assert GONGClient._can_handle_query(trange, Instrument('z'))
    assert GONGClient._can_handle_query(trange, Instrument('ct'), Physobs('INTENSITY'), Wavelength(6563*u.AA))
    assert not GONGClient._can_handle_query(trange)
    
trange = Time('2016/6/4 00:00:00', '2016/6/4 00:30:00')
@pytest.mark.online
def test_query():
    qr = GONGClient.query(trange, physobs = 'LOS_MAGNETIC_FIELD', instrument='bb')
    assert len(qr) == 3
    assert qr.time_range()[0] == '2016/06/04'
    assert qr.time_range()[1] == '2016/06/04'

#Downoads 4 fits files each of
#size 1.2MB. Total size = 4.8MB.
@pytest.mark.online
@pytest.mark.parametrize("time, physobs, instrument, wavelength",
                         [(Time('2016/6/13 03:00', '2016/6/13 04:00'), Physobs('INTENSITY'),
                           Instrument('ud'), Wavelength(6768*u.AA))])
def test_get(time, physobs, instrument, wavelength):
    qr = GONGClient.query(time, physobs, instrument, wavelength)
    res = GONGClient.get(qr)
    download_list = res.wait()
    assert len(qr) == len(download_list)

#Downloads 2 fits.gz files of size 1.2MB each.
#Total size is 2.5MB.
@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2016/6/4', '2016/6/4 00:10:00'), a.Physobs('INTENSITY'), a.Wavelength(6768*u.AA))
    assert isinstance(qr, UnifiedResponse)
    response = Fido.fetch(qr)
    assert len(response) == qr._numfile
