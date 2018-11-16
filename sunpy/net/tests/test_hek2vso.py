# -*- coding: utf-8 -*-
# Author: Michael Malocha
# e-mail: mmalocha13@gmail.com
# Version: June 11th, 2013
#

"""
This module was built to test the HEK2VSO class.
"""

__author__ = 'Michael Malocha'
__version__ = 'June 11th, 2013'

import pytest

from astropy import units as u
from astropy import table
from sunpy.time import parse_time
from sunpy.net import hek
from sunpy.net import vso
from sunpy.net import hek2vso


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


@pytest.mark.remote_data
def test_translate_results_to_query():
    """Make sure that conversion of HEK results to VSO queries is accurate"""
    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    vso_query = hek2vso.translate_results_to_query(hek_query)

    # Comparing length of two lists
    assert len(hek_query) == len(vso_query)
    # Comparing types of both queries
    # Not sure this test makes any sense now
    assert isinstance(hek_query, table.Table)
    assert isinstance(vso_query, list)


@pytest.mark.remote_data
def test_vso_attribute_parse():
    """Make sure that Parsing of VSO attributes from HEK queries is accurate"""
    h = hek.HEKClient()
    hek_query = h.search(hekTime, hekEvent)
    vso_query = hek2vso.vso_attribute_parse(hek_query[0])

    # Checking Time
    assert vso_query[0].start == parse_time(hek_query[0]['event_starttime'])
    assert vso_query[0].end == parse_time(hek_query[0]['event_endtime'])

    # Checking Observatory
    assert vso_query[1].value == hek_query[0]['obs_observatory']

    # Checking Instrument
    assert vso_query[2].value == hek_query[0]['obs_instrument']

    # Checking Wavelength
    assert vso_query[3].min == hek_query[0]['obs_meanwavel'] * u.Unit(hek_query[0]['obs_wavelunit'])
    assert vso_query[3].max == hek_query[0]['obs_meanwavel'] * u.Unit(hek_query[0]['obs_wavelunit'])
    assert vso_query[3].unit == u.Unit('Angstrom')


class TestH2VClient(object):
    """Tests the H2V class"""
    # TODO
