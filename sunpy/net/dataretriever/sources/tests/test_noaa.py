import datetime
from unittest import mock

import pytest

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.net._attrs import Instrument, Time
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.dataretriever.sources import noaa
from sunpy.tests.helpers import no_vso
from sunpy.time import parse_time
from sunpy.time.timerange import TimeRange


@pytest.fixture
def LCClient():
    return noaa.NOAAIndicesClient()


@pytest.fixture
def SRSClient():
    return noaa.SRSClient()


def mock_query_object(start_date, end_date):
    """
    Creation of a QueryResponse object, and prefill some
    downloaded data from noaa.NOAAIndicesClient().fetch(Time('20 ..)
    """
    # Create a mock QueryResponse object
    map_ = {
        'TimeRange': TimeRange(parse_time(start_date), parse_time(end_date)),
        'Time_start': parse_time(start_date),
        'Time_end': parse_time(end_date),
        'source': 'sdic',
        'instrument': 'noaa-indices',
        'physobs': 'sunspot number',
        'provider': 'swpc'
    }

    resp = QueryResponse.create(map_,
                                noaa.NOAAIndicesClient()._get_url_for_timerange(None),
                                client=noaa.NOAAIndicesClient())
    return resp


@pytest.mark.remote_data
def test_fetch_working(LCClient, tmpdir):
    """
    Tests if the online server for noaa is working.
    Uses the url : https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json
    """
    qr1 = LCClient.search(Time('2012/10/4', '2012/10/6'),
                          Instrument('noaa-indices'))

    # Mock QueryResponse object
    mock_qr = mock_query_object('2012/10/4', '2012/10/6')
    # Compare if two objects have the same attribute

    mock_qr = mock_qr.blocks[0]
    qr = qr1.blocks[0]

    assert mock_qr.source == qr.source
    assert mock_qr.provider == qr.provider
    assert mock_qr.physobs == qr.physobs
    assert mock_qr.instrument == qr.instrument
    assert mock_qr.url == qr.url
    assert mock_qr.time == qr.time

    # Assert if the timerange is same
    assert qr1.time_range() == TimeRange('2012/10/4', '2012/10/6')

    target_dir = tmpdir.mkdir("down")
    download_list = LCClient.fetch(qr1, path=target_dir)
    assert len(download_list) == len(qr1)
    assert download_list[0].split('/')[-1] == 'observed-solar-cycle-indices.json'


@pytest.mark.parametrize(
    "timerange,url_start,url_end",
    [(TimeRange('1995/06/03', '1995/06/04'),
      'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json',
      'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json'),
     (TimeRange('2008/06/01', '2008/06/02'),
      'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json',
      'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json')])
def test_get_url_for_time_range(LCClient, timerange, url_start, url_end):
    urls = LCClient._get_url_for_timerange(timerange)
    assert isinstance(urls, list)
    assert urls[0] == url_start
    assert urls[-1] == url_end


def test_can_handle_query():
    ans1 = noaa.NOAAIndicesClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('noaa-indices'))
    assert ans1
    ans2 = noaa.NOAAIndicesClient._can_handle_query(
        Time('2012/7/7', '2012/7/7'))
    assert not ans2
    ans3 = noaa.NOAAIndicesClient._can_handle_query(
        Time('2012/8/9', '2012/8/10'), Instrument('eve'))
    assert not ans3


@mock.patch('sunpy.net.dataretriever.sources.noaa.NOAAIndicesClient.search',
            return_value=mock_query_object('2012/8/9', '2012/8/10'))
def test_query(mock_search, LCClient):
    qr1 = LCClient.search(
        Time('2012/8/9', '2012/8/10'), Instrument('noaa-indices'))
    assert isinstance(qr1, QueryResponse)
    assert len(qr1) == 1
    assert qr1.time_range().start == parse_time('2012/08/09')
    assert qr1.time_range().end == parse_time('2012/08/10')


@mock.patch('sunpy.net.dataretriever.sources.noaa.NOAAIndicesClient.search',
            return_value=mock_query_object('2012/10/4', '2012/10/6'))
# The return value of download is irrelevant
@mock.patch('parfive.Downloader.download',
            return_value=None)
@mock.patch('parfive.Downloader.enqueue_file')
def test_fetch(mock_wait, mock_search, mock_enqueue, tmp_path, LCClient):
    path = tmp_path / "sub"
    path.mkdir()
    qr1 = LCClient.search(Time('2012/10/4', '2012/10/6'),
                          Instrument('noaa-indices'))
    LCClient.fetch(qr1, path=path / "{file}")

    # Here we assert that the `fetch` function has called the parfive
    # Downloader.enqueue_file method with the correct arguments. Everything
    # that happens after this point should either be tested in the
    # GenericClient tests or in parfive itself.
    assert mock_enqueue.called_once_with(("https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json",
                                          path / "observed-solar-cycle-indices.json"))


@mock.patch('sunpy.net.dataretriever.sources.noaa.NOAAIndicesClient.search',
            return_value=mock_query_object('2012/10/4', '2012/10/6'))
# The return value of download is irrelevant
@mock.patch('parfive.Downloader.download',
            return_value=None)
@mock.patch('parfive.Downloader.enqueue_file')
@no_vso
def test_fido(mock_wait, mock_search, mock_enqueue, tmp_path, LCClient):
    path = tmp_path / "sub"
    path.mkdir()
    qr1 = Fido.search(Time('2012/10/4', '2012/10/6'),
                      Instrument('noaa-indices'))
    Fido.fetch(qr1, path=path)

    # Here we assert that the `fetch` function has called the parfive
    # Downloader.enqueue_file method with the correct arguments. Everything
    # that happens after this point should either be tested in the
    # GenericClient tests or in parfive itself.
    assert mock_enqueue.called_once_with(("https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json",
                                          path / "observed-solar-cycle-indices.json"))


@no_vso
@pytest.mark.remote_data
def test_srs_tar_unpack():
    qr = Fido.search(a.Instrument("soon") & a.Time("2015/01/01", "2015/01/01T23:59:29"))
    res = Fido.fetch(qr)
    assert len(res) == 1
    assert res.data[0].endswith("20150101SRS.txt")


@no_vso
@pytest.mark.remote_data
def test_srs_tar_unpack_midyear():
    qr = Fido.search(a.Instrument("soon") & a.Time("2011/06/07", "2011/06/08T23:59:29"))
    res = Fido.fetch(qr)
    assert len(res) == 2
    files = sorted(res.data)
    assert files[0].endswith("20110607SRS.txt")
    assert files[-1].endswith("20110608SRS.txt")


@no_vso
@pytest.mark.remote_data
def test_srs_current_year():
    year = datetime.date.today().year
    qr = Fido.search(a.Instrument("soon") & a.Time(f"{year}/01/01", f"{year}/01/01T23:59:29"))
    res = Fido.fetch(qr)
    assert len(res) == 1
    assert res.data[0].endswith(f"{year}0101SRS.txt")


@no_vso
@pytest.mark.remote_data
def test_srs_save_path(tmpdir):
    qr = Fido.search(a.Instrument.srs_table, a.Time("2016/10/01", "2016/10/02"))
    files = Fido.fetch(qr, path=str(tmpdir))
    assert len(files) == 2
    files = sorted(files)
    assert files[0].endswith("20161001SRS.txt")
    assert files[1].endswith("20161002SRS.txt")


@pytest.mark.remote_data
@pytest.mark.filterwarnings('ignore:ERFA function')
def test_srs_out_of_range(SRSClient):
    res = SRSClient.search(a.Time('1995/01/01', '1995/02/01'))
    assert len(res) == 0
    res = SRSClient.search(a.Time('2995/01/01', '2995/02/01'))
    assert len(res) == 0


@pytest.mark.remote_data
@pytest.mark.filterwarnings('ignore:ERFA function')
def test_srs_start_or_end_out_of_range(SRSClient):
    res = SRSClient.search(a.Time('1995/12/30', '1996/01/02'))
    assert len(res) == 1
    cur_year = datetime.date.today().year
    res = SRSClient.search(a.Time(f'{cur_year}/01/01', f'{cur_year+2}/01/01'))
    assert len(res) > 0


@pytest.mark.remote_data
def test_tar_file_broken():
    # 2010 extracts out to 2010_SRS while other years do SRS only.
    results = Fido.search(a.Time("2010/5/1", "2010/5/2"), a.Instrument.soon)
    results = Fido.fetch(results)
    assert len(results) == 2


def test_attr_reg():
    assert a.Instrument.noaa_indices == a.Instrument("NOAA-Indices")
    assert a.Instrument.noaa_predict == a.Instrument("NOAA-Predict")
    assert a.Instrument.srs_table == a.Instrument("SRS-Table")
    assert a.Instrument.soon == a.Instrument("SOON")


def test_client_repr(LCClient):
    """
    Repr check
    """
    output = str(LCClient)
    assert output[:50] == 'sunpy.net.dataretriever.sources.noaa.NOAAIndicesCl'
