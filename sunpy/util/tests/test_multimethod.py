# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import, division, print_function

import pytest

from sunpy.util.multimethod import MultiMethod, FAIL, WARN, TypeWarning

def test_super():
    class String(str):
        pass

    mm = MultiMethod(lambda *a: a)

    @mm.add_dec(str, str)
    def foo(foo, bar):
        return 'String'

    @mm.add_dec(String, str)
    def foo(foo, bar):
        return 'Fancy', mm.super(super(String, foo), bar)

    assert mm('foo', 'bar') == 'String'
    assert mm(String('foo'), 'bar') == ('Fancy', 'String')


def test_override(recwarn):
    class String(str):
        pass

    mm = MultiMethod(lambda *a: a)

    @mm.add_dec(str, str)
    def foo(foo, bar):
        return 'String'

    pytest.raises(
        TypeError, mm.add_dec(String, str, override=FAIL), lambda x, y: None
    )

    mm.add_dec(String, str, override=WARN)(lambda x, y: None)
    w = recwarn.pop(TypeWarning)
    assert 'Definition (String, str) overrides prior definition (str, str).' in str(w.message)
