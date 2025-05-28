"""
These tests are borrowed from CPython
"""
import sys
import unittest
import collections

from sunpy.util.functools import seconddispatch


class TestSingleDispatch(unittest.TestCase):
    def test_simple_overloads(self):

        @seconddispatch
        def g(dummy, obj):
            return "base"

        def g_int(dummy, i):
            return "integer"
        g.register(int, g_int)
        assert g(None, "str") == "base"
        assert g(None, 1) == "integer"
        assert g(None, [1, 2, 3]) == "base"

    def test_mro(self):

        @seconddispatch
        def g(obj):
            return "base"

        class A:
            pass

        class C(A):
            pass

        class B(A):
            pass

        class D(C, B):
            pass

        def g_A(dummy, a):
            return "A"

        def g_B(dummy, b):
            return "B"
        g.register(A, g_A)
        g.register(B, g_B)
        assert g(None, A()) == "A"
        assert g(None, B()) == "B"
        assert g(None, C()) == "A"
        assert g(None, D()) == "B"

    def test_register_decorator(self):

        @seconddispatch
        def g(dummy, obj):
            return "base"

        @g.register(int)
        def g_int(dummy, i):
            return f"int {i}"
        assert g(None, "") == "base"
        assert g(None, 12) == "int 12"
        assert g.dispatch(int) is g_int
        assert g.dispatch(object) is g.dispatch(str)
        # Note: in the assert above this is not g.
        # @singledispatch returns the wrapper.

    def test_wrapping_attributes(self):
        @seconddispatch
        def g(dummy, obj):
            "Simple test"
            return "Test"
        assert g.__name__ == "g"
        if sys.flags.optimize < 2:
            assert g.__doc__ == "Simple test"

    def test_annotations(self):
        @seconddispatch
        def i(dummy, arg):
            return "base"

        @i.register
        def _(dummy, arg: collections.abc.Mapping):
            return "mapping"

        @i.register
        def _(dummy, arg: "collections.abc.Sequence"):
            return "sequence"
        assert i(None, None) == "base"
        assert i(None, {"a": 1}) == "mapping"
        assert i(None, [1, 2, 3]) == "sequence"
        assert i(None, (1, 2, 3)) == "sequence"
        assert i(None, "str") == "sequence"

        # Registering classes as callables doesn't work with annotations,
        # you need to pass the type explicitly.
        @i.register(str)
        class _:
            def __init__(self, dummy, arg):
                self.arg = arg

            def __eq__(self, other):
                return self.arg == other
        assert i(None, "str") == "str"
