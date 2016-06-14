import pytest

import astropy.units as u

import sunpy.net.vso.attrs as attrs
from sunpy.net.dataretriever.downloader_factory import Fido, UnifiedResponse

# Test that the correct client is called.
@pytest.mark.parametrize("query, client", [
    ((attrs.Time('2012/8/8', '2012/8/9'), attrs.Instrument('eve'),
      attrs.Level(0)), "EVEClient"),
    ((attrs.Time('2012/3/5', '2012/3/7'), attrs.Instrument('goes')),
     "GOESClient"),
    ((attrs.Time('2012/3/5', '2012/3/7'), attrs.Instrument('lyra')),
     "LYRAClient"),
    ((attrs.Time('2012/1/8', '2012/1/9'), attrs.Instrument('norh')),
     "NoRHClient"),
    ((attrs.Time('2012/1/8', '2012/3/9'), attrs.Instrument('noaa-indices')),
     "NOAAIndicesClient"),
    ((attrs.Time('2012/12/8', '2012/12/9'), attrs.Instrument('noaa-predict')),
     "NOAAPredictClient"),
])
def test_offline_client(query, client):
    unifiedresp = Fido.search(*query)
    for block in unifiedresp:
        assert block.client.__class__.__name__ == client


# Test that the correct client is called (online, so only do one for speed).
@pytest.mark.online
@pytest.mark.parametrize("query, client", [
    ((attrs.Time('2012/8/8', '2012/8/9'), attrs.Instrument('eve')),
     "VSOClient"),
])
def test_online_client(query, client):
    unifiedresp = Fido.search(*query)
    for block in unifiedresp:
        assert block.client.__class__.__name__ == client


@pytest.mark.online
@pytest.mark.parametrize(
    "time,instrument",
    [
     #(attrs.Time('2013/1/8','2013/1/9'),attrs.Instrument('eve')),
     (attrs.Time('2011/5/5', '2011/5/6'), attrs.Instrument('lyra')),
     (attrs.Time('2012/7/8', '2012/7/9'), attrs.Instrument('norh')),
     (attrs.Time('2012/11/27', '2012/11/27'), attrs.Instrument('rhessi')),
     #(attrs.Time('2012/1/8', '2012/3/9'),attrs.Instrument('noaa-indices')),
     (attrs.Time('2012/12/8', '2012/12/9'),
      attrs.Instrument('noaa-predict')),
    ])
def test_fetch(time, instrument):
    unifiedresp = Fido.search(time, instrument)
    res = Fido.fetch(unifiedresp, wait=False)
    download_list = res.wait()
    assert len(download_list) == unifiedresp.file_num


@pytest.mark.online
@pytest.mark.parametrize(
    "time1,time2,instrument",
    [
     (attrs.Time('2012/4/5', '2012/4/5'), attrs.Time(
         '2011/8/7', '2011/8/7'), attrs.Instrument('norh')),
     (attrs.Time('2012/3/3', '2012/3/3'), attrs.Time(
         '2013/1/1', '2013/1/1'), attrs.Instrument('rhessi')),
     (attrs.Time('2012/7/7', '2012/7/8'), attrs.Time(
         '2011/7/7', '2011/7/8'), attrs.Instrument('lyra')),
     #(attrs.Time('2012/3/2','2012/3/3'),attrs.Time('2012/1/27','2012/1/27'),attrs.Instrument('eve')),
    ])
def test_multiple_time(time1, time2, instrument):

    unifiedresp = Fido.search(time1 | time2, instrument)
    num_files_to_download = unifiedresp.file_num
    res = Fido.fetch(unifiedresp, wait=False)
    files_downloaded = len(res.wait())
    assert files_downloaded == num_files_to_download


@pytest.mark.online
@pytest.mark.parametrize(
    "time,instrument1,instrument2",
    [
        #(attrs.Time('2011/12/1','2011/12/4'),attrs.Instrument('eve'),attrs.Instrument('norh')),
        (attrs.Time('2011/3/5', '2011/3/5'), attrs.Instrument('lyra'),
         attrs.Instrument('rhessi')),
    ])
def test_multiple_clients(time, instrument1, instrument2):

    unifiedresp = Fido.search(time, instrument1 | instrument2)
    num_files_to_download = unifiedresp.file_num
    res = Fido.fetch(unifiedresp, wait=False)
    files_downloaded = len(res.wait())
    assert files_downloaded == num_files_to_download


@pytest.mark.online
def test_vso():
    unifiedresp = Fido.search(
        attrs.Time("2013/3/4 01:00:00", "2013/3/4 01:10:00"),
        attrs.Instrument('aia'), attrs.Wavelength(304 * u.AA, 304 * u.AA),
        attrs.Sample(600 * u.s))
    num_files_to_download = unifiedresp.file_num
    res = Fido.fetch(unifiedresp, wait=False)
    files_downloaded = len(res.wait())
    assert files_downloaded == num_files_to_download


@pytest.mark.xfail
def test_unifiedresponse_slicing():
    """
    This should pass, a fix is incoming from @Cadair
    """
    results = Fido.search(attrs.Time("2012/1/1", "2012/1/5"),
                          attrs.Instrument("lyra"))
    assert isinstance(results[0:2], UnifiedResponse)
    assert isinstance(results[0], UnifiedResponse)
