# This module was developed with funding provided by
# the Google Summer of Code 2016.
import pytest
import datetime

from astropy import units as u
from sunpy.net.vso.attrs import Time, Instrument, Physobs, Wavelength
from sunpy.net.fido_factory import UnifiedResponse
from sunpy.net import Fido
from sunpy.net import attrs as a
import sunpy.net.dataretriever.sources.vsm as vsm

VClient = vsm.VSMClient()

trange = Time('2014/6/1', '2014/6/4')


# Not downloading fits files since
# they are very huge ranging from 75MB-100MB
@pytest.mark.online
@pytest.mark.parametrize("time, instrument, wavelength, physobs",
                         [(trange, a.Instrument('vsm'), a.Wavelength(6302 * u.AA),
                           Physobs("LOS_MAGNETIC_FIELD"))])
def test_query(time, instrument, physobs, wavelength):
    qr = VClient.search(time, instrument, physobs, wavelength)
    assert len(qr) == 3


@pytest.mark.parametrize(
    "timerange, instrument, wavelength, physobs, expected",
    [(trange, a.Instrument('vsm'), a.Wavelength(6302 * u.AA), None, True),
     (trange, a.Instrument('vsm'), None, None, False),
     (trange, a.Instrument('vsm'), a.Wavelength(8542 * u.AA), None, True), (
         trange, Instrument('vsm'), a.Wavelength(1083.0 * u.nm),
         a.Physobs('EQUIVALENT_WIDTH'), True)])
def test_can_handle_query(timerange, instrument, wavelength, physobs,
                          expected):
    assert VClient._can_handle_query(timerange, instrument, wavelength,
                                     physobs) is expected


# Because otherwise the `._map` contains the latest Physobs from test_query (last search)
VClient2 = vsm.VSMClient()


@pytest.mark.online
def test_query_2():
    qr = VClient2.search(trange, Instrument('vsm'), Wavelength(6302 * u.AA))
    assert len(qr) == 9
    assert qr.time_range().start.date() == datetime.date(2014, 6, 1)
    assert qr.time_range().end.date() == datetime.date(2014, 6, 4)


# FDISK files for VECTOR_MAGNETIC_FIELD
# For wavelength 6302 Angstroms.
@pytest.mark.online
def test_fido_query():
    qr = Fido.search(
        a.Time('2015/1/3', '2015/1/6'),
        a.Instrument('vsm'),
        a.Wavelength(630.2 * u.nm), a.Physobs('VECTOR_MAGNETIC_FIELD'))
    assert isinstance(qr, UnifiedResponse)
    assert qr._numfile == 6
