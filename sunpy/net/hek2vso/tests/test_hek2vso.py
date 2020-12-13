# Author: Michael Malocha
# e-mail: mmalocha13@gmail.com
# Version: June 11th, 2013
#

"""
This module was built to test the HEK2VSO class.
"""

__author__ = "Michael Malocha"
__version__ = "June 11th, 2013"

import pytest

import astropy.units as u

from sunpy.net import attrs as a
from sunpy.net import hek, hek2vso, vso
from sunpy.time import parse_time

startTime = "2011/08/09 07:23:56"
endTime = "2011/08/09 12:40:29"
eventType = "FL"
instrument = "eit"

hekTime = a.Time(startTime, endTime)
hekEvent = a.hek.EventType(eventType)


@pytest.fixture(scope="function")
@pytest.mark.remote_data
def h2v_client():
    return hek2vso.H2VClient()


@pytest.fixture(scope="function")
@pytest.mark.remote_data
def hek_client():
    return hek.HEKClient()


@pytest.fixture
@pytest.mark.remote_data
def vso_client():
    vso.VSOClient()


@pytest.mark.remote_data
def test_translate_results_to_query(hek_client):
    """Make sure that conversion of HEK results to VSO queries is accurate"""
    h = hek_client
    hek_query = h.search(hekTime, hekEvent)
    vso_query = hek2vso.translate_results_to_query(hek_query)

    # Comparing length of two lists
    assert len(hek_query) == len(vso_query)
    # Comparing types of both queries
    # Not sure this test makes any sense now
    assert isinstance(hek_query, hek.HEKResponse)
    assert isinstance(vso_query, list)


@pytest.mark.remote_data
def test_vso_attribute_parse(hek_client):
    """Make sure that Parsing of VSO attributes from HEK queries is accurate"""
    h = hek_client
    hek_query = h.search(hekTime, hekEvent)
    vso_query = hek2vso.vso_attribute_parse(hek_query[0])

    # Checking Time
    assert vso_query[0].start == parse_time(hek_query[0]["event_starttime"])
    assert vso_query[0].end == parse_time(hek_query[0]["event_endtime"])

    # Checking Observatory
    assert vso_query[1].value == hek_query[0]["obs_observatory"]

    # Checking Instrument
    assert vso_query[2].value == hek_query[0]["obs_instrument"]

    # Checking Wavelength
    assert vso_query[3].min == hek_query[0]["obs_meanwavel"] * u.Unit(
        hek_query[0]["obs_wavelunit"]
    )
    assert vso_query[3].max == hek_query[0]["obs_meanwavel"] * u.Unit(
        hek_query[0]["obs_wavelunit"]
    )
    assert vso_query[3].unit == u.Unit("Angstrom")


@pytest.mark.remote_data
def test_members(h2v_client):
    client = h2v_client
    assert isinstance(client.hek_client, hek.HEKClient)
    assert isinstance(client.vso_client, vso.VSOClient)
    assert client.hek_results == ""
    assert client.vso_results == []
    assert client.num_of_records == 0


@pytest.mark.remote_data
def test_translate_and_query(h2v_client, hek_client):
    h = hek_client
    h2v = h2v_client
    q = h.search(a.Time(startTime, endTime), a.hek.EventType(eventType))
    h2v_q = h2v.translate_and_query(q)

    assert len(q) == len(h2v_q)
    assert isinstance(h2v_q, list)
    assert isinstance(h2v_q[0], vso.vso.QueryResponse)


@pytest.mark.remote_data
def test_full_query(h2v_client, hek_client):
    h2v = h2v_client
    h = hek_client
    h2v_q_1 = h2v.full_query(
        (a.Time(startTime, endTime), a.hek.EventType(eventType))
    )

    assert h2v.num_of_records > 1
    assert len(h2v.vso_results) > 1
    assert len(h2v.hek_results) > 1

    h2v._quick_clean()
    q = h.search(a.Time(startTime, endTime), a.hek.EventType(eventType))
    h2v_q_2 = h2v.translate_and_query(q)

    assert len(h2v_q_1) == len(h2v_q_2)
    assert len(h2v.hek_results) == len(q)
    assert h2v.hek_results[0].get_voevent() == q[0].get_voevent()

    for i in range(len(h2v_q_1)):
        assert len(h2v_q_1[i]) == len(h2v_q_2[i])
        if i != 2:
            assert h2v_q_1[i].total_size() == h2v_q_2[i].total_size()


@pytest.mark.remote_data
def test_quick_clean(h2v_client, hek_client):
    h2v = h2v_client
    h2v_q = h2v.full_query(
        (a.Time(startTime, endTime), a.hek.EventType(eventType))
    )

    assert h2v.num_of_records != 0
    assert len(h2v.vso_results) != 0

    h2v._quick_clean()

    assert h2v.vso_results == []
    assert h2v.num_of_records == 0
