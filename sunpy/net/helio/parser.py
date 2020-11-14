"""
This module is meant to parse the HELIO registry and return WSDL endpoints to
facilitate the interfacing between further modules and HELIO.
"""
import urllib
import xml.etree.ElementTree as EL
from contextlib import closing

from bs4 import BeautifulSoup

from sunpy.net.helio import registry_links as RL

__all__ = ['webservice_parser', 'endpoint_parser', 'wsdl_retriever']

# Lifespan in seconds before a link times-out
LINK_TIMEOUT = 3


def webservice_parser(service='HEC'):
    """
    Quickly parses important contents from HELIO registry.

    Uses the link contained in registry_links in with 'service' appended
    and scrapes the web-service links contained on that webpage.

    Parameters
    ----------
    service: str
        Indicates which particular HELIO service is used. Defaults to HEC.

    Returns
    -------
    links: list or NoneType
        List of urls to registries containing WSDL endpoints.

    Examples
    --------
    >>> from sunpy.net.helio import parser
    >>> parser.webservice_parser()  # doctest: +REMOTE_DATA
    ['http://msslkz.mssl.ucl.ac.uk/helio-hec/HelioService',
     'http://festung3.oats.inaf.it:8080/helio-hec/HelioService',
     'http://festung1.oats.inaf.it:8080/helio-hec/HelioService',
     'http://hec.helio-vo.eu/helio_hec/HelioService',
     'http://msslkz.mssl.ucl.ac.uk/helio-hec/HelioLongQueryService',
     'http://festung3.oats.inaf.it:8080/helio-hec/HelioLongQueryService',
     'http://festung1.oats.inaf.it:8080/helio-hec/HelioLongQueryService',
     'http://hec.helio-vo.eu/helio_hec/HelioLongQueryService']
    """
    link = RL.LINK + '/' + service.lower()
    xml = link_test(link)
    if xml is None:
        return None
    root = EL.fromstring(xml)
    links = []

    for interface in root.iter('interface'):
        service_type = interface.attrib
        key = list(service_type.keys())
        if len(key) > 0:
            value = service_type[key[0]]
            if value == 'vr:WebService':
                for url in interface.iter('accessURL'):
                    if url.text not in links:
                        links.append(url.text)
    return links


def endpoint_parser(link):
    """
    Takes a link to a list of endpoints and parses the WSDL links.

    Feeding 1 result from webservice_parser() into endpoint_parser() at a time
    will return a list of WSDL endpoints that are contained on the page from
    that link that was passed in.

    Parameters
    ----------
    link: str
        A url to a page containing links to WSDL files.

    Returns
    -------
    endpoints: list or NoneType
        A list containing all of the available WSDL endpoints from the passed
        in url.

    Examples
    --------
    >>> from sunpy.net.helio import parser
    >>> parser.endpoint_parser('http://msslkz.mssl.ucl.ac.uk/helio-hec/HelioService')  # doctest: +REMOTE_DATA
    ['http://helio.mssl.ucl.ac.uk/helio-hec/HelioService?wsdl',
    'http://helio.mssl.ucl.ac.uk/helio-hec/HelioService1_0?wsdl',
    'http://helio.mssl.ucl.ac.uk/helio-hec/HelioService1_0b?wsdl',
    'http://helio.mssl.ucl.ac.uk/helio-hec/HelioLongQueryService?wsdl',
    'http://helio.mssl.ucl.ac.uk/helio-hec/HelioLongQueryService1_0?wsdl',
    'http://helio.mssl.ucl.ac.uk/helio-hec/HelioLongQueryService1_1?wsdl',
    'http://helio.mssl.ucl.ac.uk/helio-hec/HelioLongQueryService1_0b?wsdl',
    'http://helio.mssl.ucl.ac.uk/helio-hec/HelioTavernaService?wsdl']

    """
    endpoint_page = link_test(link)
    if endpoint_page is None:
        return None
    soup = BeautifulSoup(endpoint_page, 'html.parser')
    endpoints = []
    for web_link in soup.find_all('a'):
        url = web_link.get('href')
        if url not in endpoints:
            endpoints.append(url.replace(":80", "", 1))
    return endpoints


def taverna_parser(link):
    """
    Takes a link to a list of endpoints and parses the taverna WSDL links.

    Takes a url to a page containing a list of endpoints, then passes that url
    to endpoint_parser(). Upon receiving the resulting list from the parser
    taverna_parser() goes through the list and finds all the WSDL links for
    the taverna web-service. It then returns a list containing the filtered
    links.

    Parameters
    ----------
    link: str
        A url to a page containing links to WSDL files.

    Returns
    -------
    taverna_links: list or NoneType
        A list containing WSDL links for a taverna web-service

    Examples
    --------
    >>> from sunpy.net.helio import parser
    >>> parser.taverna_parser('http://msslkz.mssl.ucl.ac.uk/helio-hec/HelioService')  # doctest: +REMOTE_DATA
    ['http://helio.mssl.ucl.ac.uk/helio-hec/HelioTavernaService?wsdl']

    """
    endpoints = endpoint_parser(link)
    taverna_links = []
    if endpoints is None:
        return None
    for web_link in endpoints:
        if 'Taverna' in web_link and web_link not in taverna_links:
            taverna_links.append(web_link)
    if len(taverna_links) == 0:
        return None
    return taverna_links


def link_test(link):
    """
    Just a quick function to test a link.

    Quickly checks to see if the URL is a valid link; if it is it returns the
    downloaded contents of that page.

    Parameters
    ----------
    link: str
        A string containing a URL

    Returns
    -------
    webpage: str or NoneType
        String containing the webresults

    Examples
    --------
    >>> from sunpy.net.helio import parser
    >>> result = parser.link_test('http://msslkz.mssl.ucl.ac.uk/helio-hec/HelioService')  # doctest: +REMOTE_DATA

    >>> print(parser.link_test('http://rrnx.invalid_url5523.com'))  # doctest: +REMOTE_DATA
        None
    """
    try:
        with closing(urllib.request.urlopen(link, timeout=LINK_TIMEOUT)) as fd:
            return fd.read()
    except (ValueError, urllib.error.URLError):
        return None


def wsdl_retriever(service='HEC'):
    """
    Retrieves a link to a taverna WSDL file

    This is essentially the master method, from it all the other functions get
    called and it essentially knits everything together. It gets a list of
    service links via webservice_parser(), then filters the results via
    taverna_parser(). Finally it tests all the returned taverna WSDL links
    and returns the first live taverna endpoint.

    Parameters
    ----------
    service: str
        Indicates which particular HELIO service is used. Defaults to HEC.

    Returns
    -------
    wsdl: str
        URL to a single live taverna endpoint

    Examples
    --------
    >>> from sunpy.net.helio import parser
    >>> parser.wsdl_retriever()  # doctest: +REMOTE_DATA
    'http://helio.mssl.ucl.ac.uk/helio-hec/HelioTavernaService?wsdl'

    Notes
    -----
    * Currently only support for HEC exists, but it was designed so that it
        could be expanded at a later date
    * There is a 3 second timeout lifespan on links, so there is potential for
        this function to take a while to return. Timeout duration can be
        controlled through the LINK_TIMEOUT value
    """
    service_links = webservice_parser(service=service)
    if service_links:
        for link in service_links:
            wsdl_links = taverna_parser(link)
            if wsdl_links:
                for end_point in wsdl_links:
                    if end_point and link_test(end_point):
                        return end_point
    raise ValueError("No online HELIO servers can be found.")
