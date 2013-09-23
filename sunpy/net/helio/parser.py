# -*- coding: utf-8 -*-
# Author:   Michael Malocha <mjm159@humboldt.edu>
# Last Edit:  September 13th, 2013
#
# This module was developed with funding from the GSOC 2013 summer of code
#

"""
This module is meant to parse the HELIO registry and return WSDL endpoints to
facilitate the interfacing between further modules and HELIO.
"""
from __future__ import absolute_import
import requests
import xml.etree.ElementTree as EL
from sunpy.net.helio import registry_links as RL
from bs4 import BeautifulSoup

__author__ = 'Michael Malocha'
__version__ = 'September 13th, 2013'


def webservice_parser(service='HEC'):
    """
    Quickly parses important contents from HELIO registry.
    """
    link = RL.LINK
    # This if block will allow future service additions.
    # turn into a dictionary
    link += service.lower()
    if link_test(link) is None:
        return None
    xml = requests.get(link).text
    root = EL.fromstring(xml)
    links = []
    # Add loop to parse xml text for table names under the 'relatedResource' tag.
    for interface in root.iter('interface'):
        service_type = interface.attrib
        key = service_type.keys()
        if len(key) > 0:
            value = service_type[key[0]]
            if value == 'vr:WebService':
                for url in interface.iter('accessURL'):
                    if url.text not in links:
                        links.append(url.text)
    return links


def endpoint_parser(link):
    """
    Takes a link to a list of endpoints and parses the WSDL links
    """
    endpoint_page = link_test(link)
    if endpoint_page is None:
        return None
    soup = BeautifulSoup(endpoint_page)
    endpoints = []
    for web_link in soup.find_all('a'):
        endpoints.append(web_link.get('href'))
    # for link in soup.find_all('interface'):
    #     for child in link.children:
    #         endpoints.append(child.string)
    return endpoints


def taverna_parser(link):
    """
    Takes a link to a list of endpoints and parses the taverna WSDL link
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
    """
    try:
        webpage = requests.get(link, timeout=3).text
        return webpage
    except requests.exceptions.Timeout:
        return None


def wsdl_retriever(service='HEC'):
    """
    uses webservices_parser & endpoint_parser to retrieve a taverna WSDL file
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