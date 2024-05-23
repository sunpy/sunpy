
import tempfile

import pytest
from hypothesis import given

import sunpy.net.dataretriever.sources.adapt as adapt
from sunpy.net import attrs as a
from sunpy.net.dataretriever.client import QueryResponse
from sunpy.net.tests.strategies import time_attr
from sunpy.time import parse_time


@pytest.fixture
def adapt_client():
    return adapt.ADAPTClient()


@given(time_attr())
def test_can_handle_query(time):
    # Hypothesis complains if we use the fixture
    adapt_client = adapt.ADAPTClient()
    ans1 = adapt_client._can_handle_query(time, a.Instrument.adapt)
    assert ans1 is True
    ans2 = adapt_client._can_handle_query(time, a.Instrument.adapt,
                                         a.adapt.ADAPTResolution('1'))
    assert ans2 is True
    ans3 = adapt_client._can_handle_query(time, a.Instrument.adapt,
                                         a.adapt.ADAPTResolution('1'),
                                         a.adapt.ADAPTHelioData('f'))
    assert ans3 is True
    ans4 = adapt_client._can_handle_query(time)
    assert ans4 is False
    ans5 = adapt_client._can_handle_query(time, a.Instrument.adapt, a.Provider.nso)
    assert ans5 is True
    ans6 = adapt_client._can_handle_query(time, a.Instrument.adapt,
                                         a.adapt.ADAPTLonType('0'))
    assert ans6 is True


def mock_query_object(adapt_client):
    """
    Creating a Query Response object and prefilling it with some information
    """
    start = '2019-05-25T02:00:00.00'
    end = '2019-05-25T02:00:59.999'
    obj = {
        'Start Time': parse_time(start),
        'End Time': parse_time(end),
        'Instrument': 'ADAPT',
        'Physobs': 'flux',
        'Source': 'GONG',
        'Provider': 'NSO',
        'url': ("https://gong.nso.edu/adapt/maps/gong/2019/adapt40311_03i012_201905250200_i00005600n0.fts.gz")
    }
    results = QueryResponse([obj], client=adapt_client)
    return results


@pytest.mark.remote_data
def test_fetch_working(adapt_client):
    """
    Tests if the online server is working.
    This also checks if the mock is working well.
    """
    start = '2019/05/25 02:00:00'
    end = '2019/05/26 02:00:59.999'
    tr = a.Time(start, end)
    qr = adapt_client.search(tr, a.Instrument.adapt)[0]
    mock_qr = mock_query_object(adapt_client)[0]

    assert mock_qr['Source'] == qr['Source']
    assert mock_qr['Provider'] == qr['Provider']
    assert mock_qr['Instrument'] == qr['Instrument']
    assert mock_qr['url'] == qr['url']
    assert qr['Start Time'].isot == mock_qr['Start Time'].isot
    assert qr['End Time'].isot == mock_qr['End Time'].isot

    with tempfile.TemporaryDirectory() as tmpdirname:
        download_list = adapt_client.fetch(qr, path=tmpdirname)
    assert len(download_list) == 1



def test_show(adapt_client):
    mock_qr = mock_query_object(adapt_client)
    qrshow0 = mock_qr.show()
    qrshow1 = mock_qr.show('Start Time', 'Instrument')
    allcols = {'Start Time', 'End Time', 'Instrument', 'Source', 'Provider', 'url'}
    assert not allcols.difference(qrshow0.colnames)
    assert qrshow1.colnames == ['Start Time', 'Instrument']
    assert qrshow0['Instrument'][0] == 'ADAPT'


def test_attr_reg():
    assert a.Instrument.adapt == a.Instrument('ADAPT')


def test_client_repr(adapt_client):
    """
    Repr check
    """
    output = str(adapt_client)
    assert output[:50] == 'sunpy.net.dataretriever.sources.adapt.ADAPTClient\n'
