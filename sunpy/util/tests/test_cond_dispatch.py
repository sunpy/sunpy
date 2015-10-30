# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import, division, print_function

import pytest

from sunpy.util.cond_dispatch import ConditionalDispatch


def pytest_funcarg__oddeven(request):
    f = ConditionalDispatch()
    # Multiply even numbers by two.
    f.add(lambda x: 2 * x, lambda x: x % 2 == 0)
    # Multiply odd numbers by three.
    f.add(lambda x: 3 * x, lambda x: x % 2 == 1)
    return f


def test_dispatch(oddeven):
    assert oddeven(2) == 4
    assert oddeven(3) == 9


def test_wrong_sig(oddeven):
    with pytest.raises(TypeError) as exc_info:
        oddeven(y=2)
    assert "There are no functions matching your input parameter signature." in str(exc_info.value)


def test_nocond():
    f = ConditionalDispatch()
    # Multiply even numbers by two.
    f.add(lambda x: 2 * x, lambda x: x % 2 == 0)
    with pytest.raises(TypeError) as exc_info:
        f(3)
        assert "Your input did not fulfill the condition for any function." in str(exc_info.value)


def test_else():
    f = ConditionalDispatch()
    # Multiply even numbers by two.
    f.add(lambda x: 2 * x, lambda x: x % 2 == 0)
    f.add(lambda x: 3 * x)
    assert f(2) == 4
    assert f(3) == 9


def test_else2():
    # This verifies else branches do not catch cases that are covered
    # by cases added later.
    f = ConditionalDispatch()
    # Because gcd(2, 3) == 1, 2 | x and 3 | x are mutually exclusive.
    f.add(lambda x: 2 * x, lambda x: x % 2 == 0)
    f.add(lambda x: 3 * x)
    f.add(lambda x: 4 * x, lambda x: x % 3 == 0)
    assert f(2) == 4
    assert f(3) == 12
    assert f(5) == 15


def test_types():
    f = ConditionalDispatch()
    f.add(lambda x: 2 * x, lambda x: x % 2 == 0, [int])
    with pytest.raises(TypeError):
        f(2.0)
