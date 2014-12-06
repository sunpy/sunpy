# -*- coding: utf-8 -*-
# Author: Michael Malocha, Rajul
# e-mail: mmalocha13@gmail.com
# Version: December 6th, 2014
#

"""
This module was built to test the HEK2VSO class.
"""

__author__ = ['Michael Malocha', 'Rajul']
__version__ = 'December 6th, 2014'

import pytest
import datetime

from astropy import units as u

from sunpy.net import hek
from sunpy.net import vso
from sunpy.net import hek2vso

from sunpy.time import parse_time


startTime = '2011/08/09 07:23:56'
endTime = '2011/08/09 12:40:29'
eventType = 'FL'
instrument = 'eit'

hekTime = hek.attrs.Time(startTime, endTime)
hekEvent = hek.attrs.EventType(eventType)

@pytest.fixture
def h2v_client():
    return hek2vso.H2VClient()

@pytest.fixture
def hek_client():
    return hek.HEKClient()

@pytest.fixture
def vso_client():
    vso.VSOClient()

@pytest.mark.online
def test_translate_results_to_query_list(hek_client):
    """Make sure that conversion of HEK results to VSO queries is accurate"""
    hek_query = hek_client.query(hekTime, hekEvent)
    vso_query = hek2vso.translate_results_to_query(hek_query)

    # Asserting length and types of two lists
    assert len(hek_query) == len(vso_query)
    assert type(hek_query) == type(vso_query)
    assert len(vso_query) == 19
    assert isinstance(vso_query, list)

    # Asserting the type and length of individual results
    assert isinstance(vso_query[0], list)
    assert len(vso_query[0]) == 4

@pytest.mark.online
def test_translate_result_to_query_single(hek_client):
    hek_query = hek_client.query(hekTime, hekEvent)[0]
    vso_query = hek2vso.translate_results_to_query(hek_query)
    
    # Asserting result to be list and have four elements
    # For Time, Source, Instrument, Wave objects
    assert isinstance(vso_query, list)
    assert isinstance(vso_query[0], list)
    assert len(vso_query) == 1
    assert len(vso_query[0]) == 4

@pytest.mark.online
def test_vso_attribute_parse(hek_client):
    """Make sure that Parsing of VSO attributes from HEK queries is accurate"""
    hek_query = hek_client.query(hekTime, hekEvent)
    vso_query = hek2vso.vso_attribute_parse(hek_query[0])

    # Cheking Time
    assert vso_query[0].start == parse_time(hek_query[0]['event_starttime'])
    assert vso_query[0].end == parse_time(hek_query[0]['event_endtime'])
    assert vso_query[0].start == datetime.datetime(2011, 8, 8, 1, 30, 4)
    assert vso_query[0].end == datetime.datetime(2011, 8, 10, 0, 0, 4)

    # Checking Observatory
    assert vso_query[1].value == hek_query[0]['obs_observatory']
    assert vso_query[1].value == 'SDO'

    # Checking Instrument
    assert vso_query[2].value == hek_query[0]['obs_instrument']
    assert vso_query[2].value == 'AIA'

    # Checking Wavelength
    assert vso_query[3].min == hek_query[0]['obs_meanwavel'] * u.Unit(hek_query[0]['obs_wavelunit'])
    assert vso_query[3].max == hek_query[0]['obs_meanwavel'] * u.Unit( hek_query[0]['obs_wavelunit'])
    assert vso_query[3].min.round() == 211.0
    assert vso_query[3].max.round() == 211.0
    assert vso_query[3].unit == u.Unit('Angstrom')

@pytest.mark.online
def test_H2VClient_instance(h2v_client):
    assert h2v_client.hek_results == ''
    assert h2v_client.vso_results == []
    assert h2v_client.num_of_records == 0
    assert isinstance(h2v_client.hek_client, hek.HEKClient)
    assert isinstance(h2v_client.vso_client, vso.VSOClient)

@pytest.mark.online
def test_H2VClient_full_query(h2v_client):
    vso_results = h2v_client.full_query((hekTime, hekEvent))

    assert len(h2v_client.vso_results) == len(h2v_client.hek_results)
    assert len(h2v_client.vso_results) == 19
    assert len(vso_results) == len(h2v_client.hek_results)
    assert len(vso_results) == 19

@pytest.mark.online
def test_H2VClient_translate_and_query(h2v_client):
    hek_results = h2v_client.hek_client.query(hekTime, hekEvent)
    vso_results = h2v_client.translate_and_query(hek_results)

    assert len(h2v_client.vso_results) == len(hek_results)
    assert len(h2v_client.vso_results) == 19
    assert len(vso_results) == len(hek_results)
    assert len(vso_results) == 19

@pytest.mark.online
def test_H2VClient_full_query_limit(h2v_client):
    vso_results = h2v_client.full_query((hekTime, hekEvent), limit=100)
    total_records = sum([len(vso_result) for vso_result in vso_results])
    total_records_n_minus_1 = total_records - len(vso_results[-1])
    expected_limit = 100

    assert len(vso_results) == 1
    assert total_records >= expected_limit 
    assert total_records_n_minus_1 <= expected_limit

@pytest.mark.online
def test_H2VClient_translate_and_query_limit(h2v_client):
    hek_results = h2v_client.hek_client.query(hekTime, hekEvent)
    vso_results = h2v_client.translate_and_query(hek_results, limit=100)
    total_records = sum([len(vso_result) for vso_result in vso_results])
    total_records_n_minus_1 = total_records - len(vso_results[-1])
    expected_limit = 100

    assert len(vso_results) == 1
    assert total_records >= expected_limit
    assert total_records_n_minus_1 <= expected_limit
