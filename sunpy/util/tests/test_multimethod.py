# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import, division, print_function

import pytest
import sys

from sunpy.util.multimethod import MultiMethod, FAIL, WARN, TypeWarning


def test_super():
    class String(str):
        pass

    mm = MultiMethod(lambda *a: a)

    with pytest.raises(TypeError):
        mm.super()

    @mm.add_dec(str, str)
    def foo(foo, bar):
        return 'String'

    # Suppress pep8 warning "F811 redefinition of unused 'foo' ..."
    @mm.add_dec(String, str)  # noqa
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

    # Illegal value for 'override'
    pytest.raises(
        ValueError, mm.add_dec(String, String, override=sys.maxsize), lambda x, y: None
    )


def test_invalid_arg_for_override_to_add_method():

    def dummy_validator():
        pass

    mm = MultiMethod(lambda *a: a)

    # See #1991
    with pytest.raises(ValueError):
        mm.add(dummy_validator, tuple(), override=sys.maxsize)


def test_call_cached():

    def sum_together(first, second):
        return first + second

    # Setup
    mm = MultiMethod(lambda *a: a)
    mm.add(sum_together, (int, int, ))
    assert mm(3, 4) == 7

    # The following call 'should' use the previously cached call
    # You'll only see the result of this in the code coverage test though
    assert mm(1, 2) == 3


def test_no_registered_methods():

    mm = MultiMethod(lambda *a: a)

    with pytest.raises(TypeError):
        mm(2)
