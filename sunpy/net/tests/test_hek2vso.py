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
