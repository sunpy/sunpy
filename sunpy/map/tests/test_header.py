# -*- coding: utf-8 -*-

import sunpy
import sunpy.map
#==============================================================================
# Test, read, get_header and write through the file independent layer
#==============================================================================
class TestMapMeta(object):

    def test_upcasing(self):
        meta = sunpy.map.MapMeta({'wibble':1, 'WOBBLE':2})
        #__getitem__
        assert meta['wibble'] == meta['WIBBLE']
        #get
        assert meta.get('wibble') == meta.get('WIBBLE')
        #has_key
        assert ('wibble' in meta) == ('WIBBLE' in meta)
        #Copy
        meta2 = meta.copy()
        assert meta2 == meta
        #pop
        assert meta.pop('wibble') == meta2.pop('WIBBLE')
        #update
        meta.update({'spam':'eggs'})
        meta2.update({'SPAM':'eggs'})
        assert meta == meta2
        #setdefault
        meta.setdefault('dave',3)
        meta2.setdefault('DAVE',3)
        assert meta.get('DAVE') == meta2.get('dave')
        #__setitem__
        meta['wibble'] = 10
        assert meta['wibble'] == 10
        meta['WIBBLE'] = 20
        assert meta['wibble'] == 20
        #__contains__
        assert 'wibble' in meta
        assert 'WIBBLE' in meta
