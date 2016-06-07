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
FClient = gong.FARSIDEClient()


#downloads 4 fits files each of size 3MB.
#Total size = 12MB.
@pytest.mark.online
@pytest.mark.parametrize("time, physobs, instrument, wavelength",
                         [(Time('2014/6/4 00:00:00', '2014/6/4 00:07:00'),
                           Physobs('intensity'), Instrument('M'), Wavelength(6563*u.AA))])
def test_query(time, physobs, instrument, wavelength):
    qr = GONGClient.query(time, physobs, instrument, wavelength)
    res = GONGClient.get(qr)
    download_list = res.wait()
    assert len(download_list) == len(qr)
                           
