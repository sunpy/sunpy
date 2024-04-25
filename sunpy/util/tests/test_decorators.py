import warnings

import pytest

import astropy.units as u

from sunpy.util.decorators import (
    active_contexts,
    deprecate_positional_args_since,
    deprecated,
    get_removal_version,
    sunpycontextmanager,
)
from sunpy.util.exceptions import SunpyDeprecationWarning, SunpyPendingDeprecationWarning


def test_removal_version_since_lts():
    major, minor = get_removal_version('2.0')
    assert major == 2
    assert minor == 1


def test_removal_version_not_since_lts():
    major, minor = get_removal_version('2.1')
    assert major == 3
    assert minor == 1


@pytest.mark.parametrize(
    ('since', 'pending', 'warning', 'message', 'warning_message'),
    [
        ('2.0', False, SunpyDeprecationWarning, '',
         'The foo function is deprecated and may be removed in version 2.1.'),
        ('2.1', False, SunpyDeprecationWarning, '',
         'The foo function is deprecated and may be removed in version 3.1.'),
        ('2.0', True, SunpyPendingDeprecationWarning, '',
         'The foo function will be deprecated in version 2.0.'),
        ('2.0', False, SunpyDeprecationWarning,
         'Custom deprecation message', 'Custom deprecation message'),
    ]
)
def test_deprecated_warning_message(since, pending, warning, message, warning_message):
    @deprecated(since, pending=pending, message=message)
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
    assert not active_contexts.get('somefunc' , False)

    with somefunc():
        # Check that the context is active while inside
        assert active_contexts.get('somefunc', False)

    # Check that the context is not active after exiting
    assert not active_contexts.get('somefunc' , False)
