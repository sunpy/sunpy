"""
These tests are borrowed from CPython
"""
import sys
import unittest
import collections

import pytest

from sunpy.util.functools import seconddispatch


class TestSingleDispatch(unittest.TestCase):
    def test_simple_overloads(self):

        @seconddispatch
        def g(dummy, obj):
            return "base"

        def g_int(dummy, i):
            return "integer"
        g.register(int, g_int)
        self.assertEqual(g(None, "str"), "base")
        self.assertEqual(g(None, 1), "integer")
        self.assertEqual(g(None, [1, 2, 3]), "base")

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
        self.assertEqual(g(None, A()), "A")
        self.assertEqual(g(None, B()), "B")
        self.assertEqual(g(None, C()), "A")
        self.assertEqual(g(None, D()), "B")

    def test_register_decorator(self):

        @seconddispatch
        def g(dummy, obj):
            return "base"

        @g.register(int)
        def g_int(dummy, i):
            return "int %s" % (i,)
        self.assertEqual(g(None, ""), "base")
        self.assertEqual(g(None, 12), "int 12")
        self.assertIs(g.dispatch(int), g_int)
        self.assertIs(g.dispatch(object), g.dispatch(str))
        # Note: in the assert above this is not g.
        # @singledispatch returns the wrapper.

    def test_wrapping_attributes(self):
        @seconddispatch
        def g(dummy, obj):
            "Simple test"
            return "Test"
        self.assertEqual(g.__name__, "g")
        if sys.flags.optimize < 2:
            self.assertEqual(g.__doc__, "Simple test")

    @pytest.mark.skipif(sys.version_info < (3, 7), reason="requires python3.7 or higher")
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
        self.assertEqual(i(None, None), "base")
        self.assertEqual(i(None, {"a": 1}), "mapping")
        self.assertEqual(i(None, [1, 2, 3]), "sequence")
        self.assertEqual(i(None, (1, 2, 3)), "sequence")
        self.assertEqual(i(None, "str"), "sequence")

        # Registering classes as callables doesn't work with annotations,
        # you need to pass the type explicitly.
        @i.register(str)
        class _:
            def __init__(self, dummy, arg):
                self.arg = arg

            def __eq__(self, other):
                return self.arg == other
        self.assertEqual(i(None, "str"), "str")
