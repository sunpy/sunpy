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

def test_wave_unit_catcher():
    """Make sure that inter-unit conversion of wavelengths is accurate"""
    # Implementing the example test cases
    assert hek2vso.wave_unit_catcher(2.11e-06, 'cm') == 210.99999999999997 
    assert hek2vso.wave_unit_catcher(9.4e-07, 'cm') == 93.99999999999999
    assert hek2vso.wave_unit_catcher(5e-08, 'mm') == 0.4999999999999999

def test_translate_results_to_query():
    h = hek.HEKClient()
    h2v = hek2vso.H2VClient()
    q = h.query(hek.attrs.Time('2011/08/09 07:23:56', '2011/08/09 12:40:29'), hek.attrs.EventType('FL'))
    ## TODO: Finish!!

class TestH2VClient:
    """Tests the H2V class"""

