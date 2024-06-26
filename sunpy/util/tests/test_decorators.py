import warnings

import pytest

import astropy.units as u

from sunpy.util.decorators import (
    ACTIVE_CONTEXTS,
    deprecate_positional_args_since,
    deprecated,
    sunpycontextmanager,
)
from sunpy.util.exceptions import SunpyDeprecationWarning


@pytest.mark.parametrize(
    ('since', 'warning', 'message', 'warning_message'),
    [
        ('2.0', SunpyDeprecationWarning, '',
         'The foo function is deprecated and may be removed in a future version'),
        ('2.1', SunpyDeprecationWarning, '',
         'The foo function is deprecated and may be removed in a future version'),
        ('2.0', SunpyDeprecationWarning,
         'Custom deprecation message', 'Custom deprecation message'),
    ]
)
def test_deprecated_warning_message(since, warning, message, warning_message):
    @deprecated(since, message=message)
    def foo():
        pass
    with pytest.warns(warning, match=warning_message):
        warnings.simplefilter('always')
        foo()

def test_deprecate_positional_args_warns_for_function_version():
    @deprecate_positional_args_since(since="0.26")
    def f1(a, *, b):
        pass

    with pytest.warns(SunpyDeprecationWarning, match=r"From version 0.26 passing these as positional"):
        f1(1, 2)

def test_deprecate_positional_args_warns_quantity_input():
    # It has to be this order otherwise, it will not work
    @deprecate_positional_args_since(since="0.26")
    @u.quantity_input
    def f1(a, *, b: u.percent = None):
        pass

    with pytest.warns(SunpyDeprecationWarning, match=r"From version 0.26 passing these as positional"):
        f1(1, 2 * u.percent)


def test_deprecate_positional_args_warns_for_class():

    class A1:
        @deprecate_positional_args_since(since="0.26")
        def __init__(self, a, b, *, c=1, d=1):
            pass

    with pytest.warns(SunpyDeprecationWarning, match=r"Pass c=3 as keyword args"):
        A1(1, 2, 3)

    with pytest.warns(SunpyDeprecationWarning, match=r"Pass c=3, d=4 as keyword args"):
        A1(1, 2, 3, 4)
    class A2:
        @deprecate_positional_args_since(since="0.26")
        def __init__(self, a=1, b=1, *, c=1, d=1):
            pass

    with pytest.warns(SunpyDeprecationWarning, match=r"Pass c=3 as keyword args"):
        A2(1, 2, 3)

    with pytest.warns(SunpyDeprecationWarning, match=r"Pass c=3, d=4 as keyword args"):
        A2(1, 2, 3, 4)


@sunpycontextmanager
def somefunc():
    print("Entering")
    yield
    print("Exiting")


def test_somefunc_context():
    # Check that the context is not active before entering
    assert not ACTIVE_CONTEXTS.get('somefunc' , False)

    with somefunc():
        # Check that the context is active while inside
        assert ACTIVE_CONTEXTS.get('somefunc', False)

    # Check that the context is not active after exiting
    assert not ACTIVE_CONTEXTS.get('somefunc' , False)
