import os
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


@pytest.fixture
def indices_client():
    return noaa.NOAAIndicesClient()


@pytest.fixture
def predict_client():
    return noaa.NOAAPredictClient()


@pytest.fixture
def srs_client():
    return noaa.SRSClient()


def mock_query_object(start_date, end_date):
    """
    Creation of a QueryResponse object, and prefill some
    downloaded data from noaa.NOAAIndicesClient().fetch(Time('20 ..)
    """
    # Create a mock Query Response object
    start = parse_time(start_date)
    end = parse_time(end_date)
    obj = {
        'Start Time': parse_time(start),
        'End Time': parse_time(end),
        'Instrument': 'NOAA-Indices',
        'Physobs': 'sunspot number',
        'Source': 'SIDC',
        'Provider': 'SWPC',
        'url': 'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json'
    }
    results = QueryResponse([obj], client=noaa.NOAAIndicesClient())
    return results


@pytest.mark.remote_data
def test_fetch_working(indices_client, tmpdir):
    """
    Tests if the online server for noaa is working.
    Uses the url : https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json
    """
    qr1 = indices_client.search(Time('2012/10/4', '2012/10/6'),
                                Instrument('noaa-indices'))
    # Mock QueryResponse object
    mock_qr = mock_query_object('2012/10/4', '2012/10/6')
    # Compare if two objects have the same attribute
    mock_qr = mock_qr[0]
    qr = qr1[0]
    assert mock_qr['Source'] == qr['Source']
    assert mock_qr['Provider'] == qr['Provider']
    assert mock_qr['Physobs'] == qr['Physobs']
    assert mock_qr['Instrument'] == qr['Instrument']
    assert mock_qr['url'] == qr['url']
    target_dir = tmpdir.mkdir("down")
    download_list = indices_client.fetch(qr1, path=target_dir)
    assert len(download_list) == len(qr1)
    assert os.path.basename(download_list[0]) == 'observed-solar-cycle-indices.json'


@pytest.mark.parametrize(
    ('timerange', 'url_start', 'url_end'),
    [(Time('1995/06/03', '1995/06/04'),
      'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json',
      'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json'),
     (Time('2008/06/01', '2008/06/02'),
      'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json',
      'https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json')])
def test_get_url_for_time_range(indices_client, timerange, url_start, url_end):
    resp = indices_client.search(timerange)
    urls = [i['url'] for i in resp]
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
def test_query(mock_search, indices_client):
    qr1 = indices_client.search(
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
def test_fetch(mock_wait, mock_search, mock_enqueue, tmp_path, indices_client):
    path = tmp_path / "sub"
    path.mkdir()
    qr1 = indices_client.search(Time('2012/10/4', '2012/10/6'),
                                Instrument('noaa-indices'))
    indices_client.fetch(qr1, path=path / "{file}")

    # Here we assert that the `fetch` function has called the parfive
    # Downloader.enqueue_file method with the correct arguments. Everything
    # that happens after this point should either be tested in the
    # GenericClient tests or in parfive itself.
    mock_enqueue.assert_called_once_with(
        Time("2012-10-04 00:00:00.000", "2012-10-06 00:00:00.000"),
        Instrument("noaa-indices")
    )


@mock.patch('sunpy.net.dataretriever.sources.noaa.NOAAIndicesClient.search',
            return_value=mock_query_object('2012/10/4', '2012/10/6'))
# The return value of download is irrelevant
@mock.patch('parfive.Downloader.download',
            return_value=None)
@mock.patch('parfive.Downloader.enqueue_file')
@no_vso
def test_fido(mock_wait, mock_search, mock_enqueue, tmp_path, indices_client):
    path = tmp_path / "sub"
    path.mkdir()
    qr1 = Fido.search(Time('2012/10/4', '2012/10/6'),
                      Instrument('noaa-indices'))
    Fido.fetch(qr1, path=path)

    # Here we assert that the `fetch` function has called the parfive
    # Downloader.enqueue_file method with the correct arguments. Everything
    # that happens after this point should either be tested in the
    # GenericClient tests or in parfive itself.
    mock_enqueue.assert_called_once_with(
        Time('2012-10-04 00:00:00.000', '2012-10-06 00:00:00.000'),
        Instrument("noaa-indices")
    )


@no_vso
@pytest.mark.remote_data
def test_srs_unpack():
    qr = Fido.search(a.Instrument("soon") & a.Time("2015/01/01", "2015/01/01T23:59:29"))
    res = Fido.fetch(qr)
    assert len(res) == 1
    assert res.data[0].endswith("20150101SRS.txt")


@no_vso
@pytest.mark.remote_data
def test_srs_midyear():
    qr = Fido.search(a.Instrument("soon") & a.Time("2011/06/07", "2011/06/08T23:59:29"))
    res = sorted(Fido.fetch(qr))
    assert len(res) == 2
    assert res[0].endswith("20110607SRS.txt")
    assert res[-1].endswith("20110608SRS.txt")


@no_vso
@pytest.mark.remote_data
def test_srs_current_year():
    # Current year is nothing but text files, all older years should be tar files.
    year = datetime.date.today().year
    qr = Fido.search(a.Instrument("soon") & a.Time(f"{year}/01/01", f"{year}/01/01T23:59:29"))
    res = Fido.fetch(qr)
    assert len(res) <= 1
    if len(res):
        assert res.data[0].endswith(f"{year}0101SRS.txt")


@no_vso
@pytest.mark.remote_data
def test_srs_save_path(tmpdir):
    qr = Fido.search(a.Instrument.srs_table, a.Time("2016/10/01", "2016/10/02"))
    files = sorted(Fido.fetch(qr, path=str(tmpdir)))
    assert len(files) == 2
    assert files[0].endswith("20161001SRS.txt")
    assert files[1].endswith("20161002SRS.txt")


@pytest.mark.remote_data
@pytest.mark.filterwarnings('ignore:ERFA function')
def test_srs_out_of_range(srs_client):
    res = srs_client.search(a.Time('1995/01/01', '1995/02/01'))
    assert len(res) == 0
    res = srs_client.search(a.Time('2995/01/01', '2995/02/01'))
    assert len(res) == 0


@pytest.mark.remote_data
@pytest.mark.filterwarnings('ignore:ERFA function')
def test_srs_start_or_end_out_of_range(srs_client):
    res = srs_client.search(a.Time('1995/12/30', '1996/01/02'))
    assert len(res) == 1
    cur_year = datetime.date.today().year
    # Will fail on the first day of the next year.
    res = srs_client.search(a.Time(f'{cur_year}/01/01', f'{cur_year+2}/01/01'))
    assert len(res) > 0


@pytest.mark.remote_data
def test_tar_file_broken():
    # 2010 extracts out to 2010_SRS while other years do SRS only.
    results = Fido.search(a.Time("2010/5/1", "2010/5/2"), a.Instrument.soon)
    results = Fido.fetch(results)
    assert len(results) == 2


def test_no_time(predict_client, indices_client):
    res = indices_client.search(a.Instrument.noaa_indices)
    assert len(res) == 1
    res = predict_client.search(a.Instrument.noaa_predict)
    assert len(res) == 1


def test_attr_reg():
    assert a.Instrument.noaa_indices == a.Instrument("NOAA-Indices")
    assert a.Instrument.noaa_predict == a.Instrument("NOAA-Predict")
    assert a.Instrument.srs_table == a.Instrument("SRS-Table")
    assert a.Instrument.soon == a.Instrument("SOON")


def test_client_repr(indices_client):
    """
    Repr check
    """
    output = str(indices_client)
    assert output[:50] == 'sunpy.net.dataretriever.sources.noaa.NOAAIndicesCl'


def test_show():
    mock_qr = mock_query_object('2012/10/4', '2012/10/6')
    qrshow0 = mock_qr.show()
    qrshow1 = mock_qr.show('Source', 'Instrument')
    allcols = {'Start Time', 'End Time', 'Instrument', 'Physobs', 'Source', 'Provider', 'url'}
    assert not allcols.difference(qrshow0.colnames)
    assert qrshow1.colnames == ['Source', 'Instrument']
    assert qrshow0['Instrument'][0] == 'NOAA-Indices'
