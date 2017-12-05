from __future__ import absolute_import

import pytest
import mock

from sunpy.net.helio import hec
from sunpy.net.helio.parser import (endpoint_parser, link_test, taverna_parser, webservice_parser,
                                    wsdl_retriever)
from sunpy.extern.six.moves import urllib


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


def test_suds_unwrapper():
    suds_output = """<?xml version="1.0" encoding="UTF-8"?>
    <S:Envelope ..... >
       <S:Body>
          <helio:queryResponse ... >
             <VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.1" version="1.1">
                <RESOURCE>
                ...
                </RESOURCE>
             </VOTABLE>
          </helio:queryResponse>
       </S:Body>
    </S:Envelope>
    """
    expected_output = """<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.1" version="1.1">
                <RESOURCE>
                ...
                </RESOURCE>
             </VOTABLE>
"""
    assert hec.suds_unwrapper(suds_output) == expected_output


@pytest.mark.remote_data
def test_webservice_parser():
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
    return ('http://helio.mssl.ucl.ac.uk:80/helio-hec/HelioTavernaService?wsdl',
            'http://helio.mssl.ucl.ac.uk:80/helio-hec/HelioLongQueryService?wsdl',
            'http://helio.mssl.ucl.ac.uk:80/helio-hec/HelioLongQueryService1_1?wsdl',
            'http://helio.ucl.ac.uk:80/helio-hec/HelioLongQueryService1_0b?wsdl')

# Test `sunpy.net.helio.parser.webservice_parser(...)`


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

# Test `sunpy.net.helio.parser.endpoint_parser(...)`


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

# `sunpy.net.helio.parser.taverna_parser(...)`


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

# Test `sunpy.net.helio.parser.wsdl_retriever(...)`


@mock.patch('sunpy.net.helio.parser.webservice_parser', return_value=None)
def test_wsdl_retriever_no_content(mock_endpoint_parser):
    """
    No links found? Return None
    """
    assert wsdl_retriever() is None


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
    Unable to find any valid Taverna URLs? Return None
    """
    assert wsdl_retriever() is None

# Test `sunpy.net.helio.parser.link_test(...)`


@mock.patch('sunpy.net.helio.parser.urllib.request.urlopen')
def test_link_test(mock_urlopen):
    """
    Read from an open, 'mocked', URL.
    """

    class MockFile(object):

        def __init__(self, content):
            self.content = content

        def read(self):
            return self.content

        def close(self):
            return

    expected = '<!doctype html><title>T</title>'
    mock_urlopen.return_value = MockFile(expected)
    assert link_test('http://www/google.com') == expected

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
