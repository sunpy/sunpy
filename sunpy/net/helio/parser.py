# -*- coding: utf-8 -*-
# Author:   Michael Malocha <mjm159@humboldt.edu>
# Last Edit:  September 7th, 2013
#
# This module was developed with funding from the GSOC 2013 summer of code
#

from __future__ import absolute_import
import requests
import xml.etree.ElementTree as EL
from sunpy.net import helio as RL
from bs4 import BeautifulSoup


def webservice_parser(service='HEC'):
    """
    Quickly parses important contents from HELIO registry
    """
    if service == 'HEC':
        link = RL.HEC
    elif service == 'HFC':
        link = RL.HFC
    else:
        print 'Unknown service'
        return None
    xml = requests.get(link).text
    root = EL.fromstring(xml)
    links = []
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
    return endpoints
    # try:
    #     endpoint_page = requests.get(link, timeout=0.5)
    #     soup = BeautifulSoup(endpoint_page)
    #     endpoints = []
    #     for link in soup.find_all('a'):
    #         endpoints.append(link.get('href'))
    #     return endpoints
    # except requests.exceptions.Timeout:
    #     return None


def link_test(link):
    """
    Just a quick function to test a link.
    """
    try:
        webpage = requests.get(link, timeout=0.5)
        return webpage
    except requests.exceptions.Timeout:
        return None


def wsdl_retriever(service='HEC'):
    """
    uses webservices_parser and endpoint_parser to retrieve a WSDL file
    """
    service_links = webservice_parser(service=service)
    if service_links is None:
        return None
    for link in service_links:
        wsdl = endpoint_parser(link)
        if wsdl is not None:
            break
    wsdl_link = None
    for link in wsdl:
        if link_test(link) is not None:
            wsdl_link = link
            break
    return wsdl_link


def main():
    """
    For some testing
    """
    print wsdl_retriever()

if __name__ == 'main':
    main()