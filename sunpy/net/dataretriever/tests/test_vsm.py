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
import sunpy.net.dataretriever.sources.vsm as vsm

VClient = vsm.VSMClient()

trange = Time('2014/6/1', '2014/6/4')
#Not downloading fits files since
#they are very huge ranging from 75MB-100MB
@pytest.mark.online
@pytest.mark.parametrize("time, instrument, wavelength, physobs",
[(trange, Instrument('vsm'), Wavelength(6302*u.AA),
Physobs("LOS_MAGNETIC_FIELD"))])
def test_query(time, instrument, physobs, wavelength):
    qr = VClient.query(time, instrument, physobs, wavelength)
    assert len(qr) == 3

def test_can_handle_query():
    assert VClient._can_handle_query(trange, Instrument('vsm'), Wavelength(6302*u.AA))
    assert not VClient._can_handle_query(trange, Instrument('vsm'))
    assert VClient._can_handle_query(trange, Instrument('vsm'), Wavelength(8542*u.AA))
    assert VClient._can_handle_query(trange, Instrument('vsm'), Wavelength(1083.0*u.nm),
                                     Physobs('EQUIVALENT_WIDTH'))

@pytest.mark.online
def test_query():
    qr = VClient.query(trange, Instrument('vsm'), Wavelength(6302*u.AA))
    assert len(qr) == 9
    assert qr.time_range()[0] == '2014/06/01'
    assert qr.time_range()[1] == '2014/06/04'

#FDISK files for VECTOR_MAGNETIC_FIELD
#For wavelength 6302 Angstroms.
@pytest.mark.online
def test_fido_query():
    qr = Fido.search(a.Time('2015/1/3', '2015/1/6'), a.Instrument('vsm'), a.Wavelength(630.2*u.nm),
                     a.Physobs('VECTOR_MAGNETIC_FIELD'))
    assert isinstance(qr, UnifiedResponse)
    assert qr._numfile == 6
