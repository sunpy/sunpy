# -*- coding: utf-8 -*-
# Author:   Michael Malocha <mjm159@humboldt.edu>
# Last Edit:  September 22nd, 2013
#
# This module was developed with funding from the GSOC 2013 summer of code
#

"""
This module is meant to parse the HELIO registry and return WSDL endpoints to
facilitate the interfacing between further modules and HELIO.
"""
from __future__ import absolute_import, print_function

import xml.etree.ElementTree as EL
from bs4 import BeautifulSoup
from contextlib import closing

from sunpy.net.helio import registry_links as RL
from sunpy.extern.six.moves import urllib

__author__ = 'Michael Malocha'
__version__ = 'September 22nd, 2013'

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
    >>> parser.webservice_parser()
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
        return xml
    root = EL.fromstring(xml)
    links = []

    #WARNING: getiterator is deprecated in Python 2.7+
    #Fix for 3.x support
    for interface in root.getiterator('interface'):
        service_type = interface.attrib
        key = list(service_type.keys())
        if len(key) > 0:
            value = service_type[key[0]]
            if value == 'vr:WebService':
                for url in interface.getiterator('accessURL'):
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
    >>> parser.endpoint_parser('http://msslkz.mssl.ucl.ac.uk/helio-hec/HelioService')
    ['http://msslkz.mssl.ucl.ac.uk:80/helio-hec/HelioService?wsdl',
    'http://msslkz.mssl.ucl.ac.uk:80/helio-hec/HelioService1_0?wsdl',
    'http://msslkz.mssl.ucl.ac.uk:80/helio-hec/HelioService1_0b?wsdl',
    'http://msslkz.mssl.ucl.ac.uk:80/helio-hec/HelioLongQueryService?wsdl',
    'http://msslkz.mssl.ucl.ac.uk:80/helio-hec/HelioLongQueryService1_0?wsdl',
    'http://msslkz.mssl.ucl.ac.uk:80/helio-hec/HelioLongQueryService1_1?wsdl',
    'http://msslkz.mssl.ucl.ac.uk:80/helio-hec/HelioLongQueryService1_0b?wsdl',
    'http://msslkz.mssl.ucl.ac.uk:80/helio-hec/HelioTavernaService?wsdl']
    """
    endpoint_page = link_test(link)
    if endpoint_page is None:
        return None
    soup = BeautifulSoup(endpoint_page)
    endpoints = []
    for web_link in soup.find_all('a'):
        endpoints.append(web_link.get('href'))
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
    >>> parser.taverna_parser('http://msslkz.mssl.ucl.ac.uk/helio-hec/HelioService')
    ['http://msslkz.mssl.ucl.ac.uk:80/helio-hec/HelioTavernaService?wsdl']
    """
    endpoints = endpoint_parser(link)
    taverna_links = []
    if endpoints is None:
        return None
    for web_link in endpoints:
        if 'Taverna' in web_link:
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
    >>> parser.link_test('http://msslkz.mssl.ucl.ac.uk/helio-hec/HelioService')
    u'<html>\n<head>...</body>\n</html>\n'

    >>> print(parser.link_test('http://rrnx.invalid_url5523.com'))
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
    >>> parser.wsdl_retriever()
    'http://msslkz.mssl.ucl.ac.uk:80/helio_hec/HelioTavernaService?wsdl'

    Notes
    -----
    * Currently only support for HEC exists, but it was designed so that it
        could be expanded at a later date
    * There is a 3 second timeout lifespan on links, so there is potential for
        this function to take a while to return. Timeout duration can be
        controlled through the LINK_TIMEOUT value
    """
    service_links = webservice_parser(service=service)
    wsdl = None
    wsdl_links = None
    if service_links is None:
        return None
    for link in service_links:
        wsdl_links = taverna_parser(link)
    if wsdl_links is None:
        return None
    for end_point in wsdl_links:
        if end_point is not None and link_test(end_point) is not None:
            wsdl = end_point
            break
    return wsdl
