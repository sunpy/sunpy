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

import unittest
from sunpy.net import hek
from sunpy.net import vso
from sunpy.net import hek2vso

startTime = '2011/08/09 07:23:56'
endTime = '2011/08/09 12:40:29'
eventType = 'FL'
instrument = 'eit'


class TestH2vClient(unittest.TestCase):
    """
    The class that runs the big happy tests
    """
    def setUp(self):
        self.h2v = hek2vso.H2VClient()
        self.hek = hek.HEKClient()
        self.vso = vso.VSOClient()
        self.hekTime = hek.attrs.Time(startTime, endTime)
        self.hekEvent = hek.attrs.EventType(eventType)

    def tearDown(self):
        self.h2v.dispose()
        self.hek.dispose()
        self.vso.dispose()
        self.hekTime.dispose()
        self.hekEvent.dispose()
        self.h2v = None
        self.hek = None
        self.vso = None
        self.hekTime = None
        self.hekEvent = None

    def test_hek_client(self):
        results = self.hek.query(self.hekTime, self.hekEvent)
        assert len(results) == 19

    def test_vso_client(self):
        result1 = self.vso.query(vso.attrs.Time(startTime, endTime),
                                 vso.attrs.Instrument(instrument))
        result2 = self.vso.query(vso.attrs.Time(startTime, endTime))
        result3 = self.vso.query(vso.attrs.Time((2011,9,20,1), (2011,9,20,2)),
                          vso.attrs.Instrument('eit'))
        self.assertEqual(len(result1), 0, 'Error in VSO query')
        self.assertEqual(len(result2), 3450, 'Error in VSO query')
        self.assertEqual(len(result3), 4, 'Error in VSO query')

    def test_translate_and_query(self):
        res = self.hek.query(hek.attrs.Time(startTime, endTime),
                             hek.attrs.EventType(eventType))
        res = self.h2v.translate_and_query(res)
        self.assertEqual(res, 'pass')

if __name__ == '__main__':
    unittest.main()