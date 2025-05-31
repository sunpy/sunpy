import urllib
from unittest import mock

import pytest
from requests.exceptions import SSLError

from sunpy.net import attrs as a
from sunpy.net.helio.hec import HECClient
from sunpy.net.helio.parser import (
    endpoint_parser,
    link_test,
    taverna_parser,
    webservice_parser,
    wsdl_retriever,
)
from sunpy.util.exceptions import SunpyUserWarning

# Currently helio makes unverified requests - this filter should be removed when
# https://github.com/sunpy/sunpy/issues/4401 is fixed
pytestmark = [pytest.mark.filterwarnings('ignore:Unverified HTTPS request is being made')]


@pytest.fixture(scope="session")
def client():
    try:
        client = HECClient()
        return client
    # If no links are found, the client should raise a ValueError
    except ValueError:
        pytest.xfail("No HELIO working links found.")


def wsdl_endpoints():
    """
    Slightly simplified form of the content on http://msslkz.mssl.ucl.ac.uk/helio-hec/HelioService
    Intentionally contains duplicate URLs
    """
    return '''
    <html><body><table>
    <tr><td>Port Name:</td><td>{http://helio-vo.eu/xml/QueryService/v1.0}HelioQueryServicePort</td>
    </tr>
    <tr><td><a href="http://helio.org/hec/HS1_0?wsdl">http://helio.org/hec/HS1_0?wsdl</a></td></tr>
    <tr><td><a href="http://helio.org/hec/HS1_0b?wsdl">http://helio.org/hec/HS1_0b?wsdl</a></td>
    </tr>
    <tr><td><a href="http://helio.org/hec/HLQS?wsdl">http://helio.org/hec/HLQS?wsdl</a></td></tr>
    <tr><td><a href="http://helio.org/hec/HLQS1_0?wsdl">http://helio.org/hec/HLQS1_0?wsdl</a></td>
    </tr>
    <tr><td><a href="http://helio.org/hec/HS1_0?wsdl">http://helio.org/hec/HS1_0?wsdl</a></td></tr>
    </table></body></html>
    '''


def hec_urls():
    """
    intentionally contains duplicate 'accessURL' elements
    """
    return '''
    <ri:Resource xmlns:ri="http://www.ivoa.net/xml/RegistryInterface/v1.0"
                 xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                 xsi:type="vs:CatalogService">
      <capability standardID="ivo://helio-vo.eu/std/FullQuery/v0.2">
        <interface xsi:type="vs:ParamHTTP">
          <accessURL use="full">http://helio.uk/hec/HelioQueryService</accessURL>
        </interface>
        <interface xsi:type="vr:WebService">
          <accessURL use="full">http://helio.uk/hec/HelioService</accessURL>
          <accessURL use="full">http://msslkk.uk/hec/HelioService</accessURL>
          <accessURL use="full">http://voparis.fr/hec/helio-hec/HelioService</accessURL>
          <accessURL use="full">http://hec.eu/helio_hec/HelioService</accessURL>
        </interface>
      </capability>
      <capability standardID="ivo://helio-vo.eu/std/FullQuery/Soap/v1.0">
        <interface xsi:type="vr:WebService">
          <accessURL use="full">http://helio.uk/hec/HelioService</accessURL>
        </interface>
      </capability>
      <capability standardID="ivo://helio-vo.eu/std/LongFullQuery/Soap/v1.0">
        <interface xsi:type="vr:WebService">
          <accessURL use="full">http://helio.uk/hec/HelioLongQueryService</accessURL>
          <accessURL use="full">http://hec.eu/helio_hec/HelioLongQueryService</accessURL>
        </interface>
      </capability>
    </ri:Resource>
    '''


@pytest.mark.remote_data
def test_webservice_parser(client):  # NOQA: ARG001
    # The client is used to check if HELIO is working
    result = webservice_parser()
    assert isinstance(result, list)


def some_taverna_urls():
    """
    Some valid `Taverna` links, duplicates intentional
    """
    return ('http://www.helio.uk/Taverna/hec?wsdl',
            'http://not.a.taverna.link/helio?wsdl',
            'http://www.abc.ord/HelioTavernaService?wsdl',
            'http://another.not.a.taverna.link/helio?wsdl',
            'http://www.helio.uk/Taverna/hec?wsdl')


def wsdl_urls():
    """
    No `Taverna` links, just `WSDL`
    """
    return ('http://helio.mssl.ucl.ac.uk/helio-hec/HelioTavernaService?wsdl',
            'http://helio.mssl.ucl.ac.uk/helio-hec/HelioLongQueryService?wsdl',
            'http://helio.mssl.ucl.ac.uk/helio-hec/HelioLongQueryService1_1?wsdl',
            'http://helio.ucl.ac.uk/helio-hec/HelioLongQueryService1_0b?wsdl')


@mock.patch('sunpy.net.helio.parser.link_test', return_value=None)
def test_webservice_parser_no_content(mock_link_test):
    """
    No content from supplied URL? Return None
    """
    assert webservice_parser('http://www.google.com') is None


@mock.patch('sunpy.net.helio.parser.link_test', return_value=hec_urls())
def test_webservice_parser_get_links(mock_link_test):
    """
    The `sunpy.net.helio.parser.link_test` returns an XML fragment with
    embedded `accessURL` elements. Ensure that all the `accessURL` are
    extracted and duplicates discarded.
    """
    hec_links = webservice_parser('http://www.google.com')

    assert len(hec_links) == 6

    assert 'http://helio.uk/hec/HelioService' in hec_links
    assert 'http://msslkk.uk/hec/HelioService' in hec_links
    assert 'http://voparis.fr/hec/helio-hec/HelioService' in hec_links
    assert 'http://hec.eu/helio_hec/HelioService' in hec_links
    assert 'http://helio.uk/hec/HelioLongQueryService' in hec_links
    assert 'http://hec.eu/helio_hec/HelioLongQueryService' in hec_links


@mock.patch('sunpy.net.helio.parser.link_test', return_value=None)
def test_endpoint_parser_no_content(mock_link_test):
    """
    No content from the supplied URL? Return None
    """
    assert endpoint_parser('http://example.com') is None


@mock.patch('sunpy.net.helio.parser.link_test', return_value=wsdl_endpoints())
def test_endpoint_parser_get_links(mock_link_test):
    """
    Get all the WSDL endpoints listed on the page of the supplied URL.
    Ensure duplicates are removed.
    """
    endpoints = endpoint_parser('http://www.google.com')

    assert len(endpoints) == 4
    assert 'http://helio.org/hec/HS1_0?wsdl' in endpoints
    assert 'http://helio.org/hec/HS1_0b?wsdl' in endpoints
    assert 'http://helio.org/hec/HLQS?wsdl' in endpoints
    assert 'http://helio.org/hec/HLQS1_0?wsdl' in endpoints


@mock.patch('sunpy.net.helio.parser.endpoint_parser', return_value=None)
def test_taverna_parser_no_content(mock_endpoint_parser):
    """
    No links at all? Return None
    """
    assert taverna_parser('http://example.com') is None


@mock.patch('sunpy.net.helio.parser.endpoint_parser', return_value=['http://try.the.pub/hec'])
def test_taverna_parser_no_taverna_links(mock_endpoint_parser):
    """
    There are some URLs but none of them Taverna URLs. Return `None`
    """
    assert taverna_parser('http://www.google.com') is None


@mock.patch('sunpy.net.helio.parser.endpoint_parser', return_value=some_taverna_urls())
def test_taverna_parser_get_taverna_links(mock_endpoint_parser):
    """
    Retrieve all the Taverna URLs
    """
    taverna_links = taverna_parser('http://www.google.com')
    assert len(taverna_links) == 2

    assert 'http://www.helio.uk/Taverna/hec?wsdl' in taverna_links
    assert 'http://www.abc.ord/HelioTavernaService?wsdl' in taverna_links


@mock.patch('sunpy.net.helio.parser.webservice_parser', return_value=None)
def test_wsdl_retriever_no_content(mock_endpoint_parser):
    """
    No links found? Raise ValueError
    """
    with pytest.raises(ValueError, match="No online HELIO servers can be found."):
        wsdl_retriever()


@mock.patch('sunpy.net.helio.parser.webservice_parser', return_value=wsdl_urls())
@mock.patch('sunpy.net.helio.parser.taverna_parser', return_value=some_taverna_urls())
@mock.patch('sunpy.net.helio.parser.link_test', return_value='some text read')
def test_wsdl_retriever_get_link(mock_link_test, mock_taverna_parser, mock_webservice_parser):
    """
    Get a Taverna link
    """
    assert wsdl_retriever() == 'http://www.helio.uk/Taverna/hec?wsdl'


@mock.patch('sunpy.net.helio.parser.webservice_parser', return_value=wsdl_urls())
@mock.patch('sunpy.net.helio.parser.taverna_parser', return_value=None)
def test_wsdl_retriever_no_taverna_urls(mock_taverna_parser, mock_webservice_parser):
    """
    Unable to find any valid Taverna URLs? Raise ValueError
    """
    with pytest.raises(ValueError, match="No online HELIO servers can be found."):
        wsdl_retriever()


@mock.patch('sunpy.net.helio.parser.link_test', return_value=None)
@mock.patch('sunpy.net.helio.parser.webservice_parser', return_value=wsdl_urls())
@mock.patch('sunpy.net.helio.parser.taverna_parser', return_value=some_taverna_urls())
def test_wsdl_retriever_wsdl(mock_taverna_parser, mock_webservice_parser, mock_link_test):
    """
    Unable to find any valid Taverna URLs? Raise ValueError
    """
    with pytest.raises(ValueError, match="No online HELIO servers can be found."):
        wsdl_retriever()


@pytest.mark.remote_data
def test_link_test():
    assert b"# SunPy Sample Data" in link_test('http://data.sunpy.org/sunpy/README.md')

# The following two tests for `link_test` have empty URLs as arguments. This is because
# when running the tests under Py2.7, I was getting the following error:
#
# "An attempt was made to connect to the internet by a test that was not marked `remote_data`"
#
# The empty URLs in no way invalidate the tests.


@mock.patch('sunpy.net.helio.parser.link_test', side_effect=ValueError)
def test_link_test_on_valueerror(mock_link_test):
    """
    If `link_test` internally raises `ValueError`, ensure it
    returns `None`
    """
    link_test('') is None


@mock.patch('sunpy.net.helio.parser.link_test', side_effect=urllib.error.URLError)
def test_link_test_on_urlerror(mock_link_test):
    """
    If `link_test` internally raises `URLError`, ensure it
    returns `None`
    """
    link_test('') is None


@mock.patch('sunpy.net.helio.parser.webservice_parser', return_value=wsdl_urls())
@mock.patch('sunpy.net.helio.parser.taverna_parser', return_value=some_taverna_urls())
@mock.patch('sunpy.net.helio.parser.link_test', return_value='some text read')
@mock.patch('sunpy.net.helio.hec.Client', side_effect=SSLError('SSL error'))
def test_ssl_verify_error(mock_webservice, mock_taverna, mock_link, mock_zeep, caplog):
    client = HECClient()
    query = client.search(a.Time('2023/02/03', '2023/02/03'))
    assert len(query) == 0
    assert "Set the 'NO_VERIFY_HELIO_SSL' environment variable disable SSL verification for Helio." in caplog.text


@pytest.mark.remote_data
def test_get_table_names(client):
    tables = client.get_table_names()
    assert len(tables) == 126
    table = tables[0][0]
    assert isinstance(table, str)
    assert table == 'timed_see_flare'


@pytest.mark.remote_data
def test_select_table(client, monkeypatch):
    monkeypatch.setattr('builtins.input', lambda x: "11")
    assert isinstance(client.select_table(), str)
    monkeypatch.setattr('builtins.input', lambda x: "e")
    assert client.select_table() is None


@pytest.mark.xfail(reason="This test is failing because the server is returning a 500 error.")
@pytest.mark.remote_data
def test_client_search(client):
    start = '2005/01/03'
    end = '2005/12/03'
    table_name = 'rhessi_hxr_flare'
    with pytest.warns(SunpyUserWarning, match="Number of results is the same as current limit. "):
        res = client.search(a.Time(start, end), a.helio.TableName(table_name), a.helio.MaxRecords(10))
    assert len(res) == 10


def test_max_records_limit():
    with pytest.raises(ValueError, match="Helio will only return a max of 20000 results."):
        a.helio.MaxRecords(99999)


@pytest.mark.xfail(reason="This test is failing because the server is returning a 500 error.")
@pytest.mark.remote_data
def test_HECResponse_iter(client):
    start = '2005/01/03'
    end = '2005/12/03'
    table_name = 'rhessi_hxr_flare'
    res = client.search(a.Time(start, end), a.helio.TableName(table_name), a.helio.MaxRecords(10000))
    for i in res:
        # Just to make sure iter still works, check number of columns
        assert len(i) == 13
