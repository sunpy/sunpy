import warnings

import pytest

from sunpy.util.decorators import _active_contexts, cached_property_based_on, deprecated, sunpycontextmanager
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
    with pytest.warns(warning, match=warning_message):  # NOQA: PT031
        warnings.simplefilter('always')
        foo()


@sunpycontextmanager
def ctx1():
    yield


@sunpycontextmanager
def ctx2():
    yield


def test_context_tracking():
    ctx1_name = f"{ctx1.__module__}.{ctx1.__qualname__}"
    ctx2_name = f"{ctx2.__module__}.{ctx2.__qualname__}"

    # Check that no sunpy contexts are active before entering
    assert _active_contexts.get() == []

    with ctx1():
        # Check that the context is active while inside
        assert _active_contexts.get() == [ctx1_name]

        with ctx2():
            # Check nesting of contexts
            assert _active_contexts.get() == [ctx1_name, ctx2_name]

            with ctx1():
                # Check a repeated context in the nesting
                assert _active_contexts.get() == [ctx1_name, ctx2_name, ctx1_name]

            # Check that only the last context is removed and not its duplicate
            assert _active_contexts.get() == [ctx1_name, ctx2_name]

        assert _active_contexts.get() == [ctx1_name]

    assert _active_contexts.get() == []


def test_cached_property_based_on():
    class Foo:
        def __init__(self, attr):
            self._attr = attr
            self._value = attr
            self.n_calls = 0

        @property
        def attr(self):
            return self._attr

        @property
        @cached_property_based_on('attr')
        def prop(self):
            self.n_calls += 1
            return self._value

    foo = Foo(1)
    assert foo.prop == 1
    assert foo.n_calls == 1

    # Changing `_value` while `attr` stays the same should not cause the
    # property to be recomputed, so `prop` must not move even though the
    # underlying value did.
    foo._value = 99
    assert foo.prop == 1
    assert foo.n_calls == 1

    # Changing `attr` should cause the property to be recomputed.
    foo._attr = 2
    foo._value = 2
    assert foo.prop == 2
    assert foo.n_calls == 2


def test_shared_attr_cache_across_properties():
    # 2 properties that both cache off the same attr must invalidate separately
    # otherwoise, reading one property would mark the other as valid incorrectly.
    class Foo:
        def __init__(self, attr):
            self._attr = attr

        @property
        def attr(self):
            return self._attr

        @property
        @cached_property_based_on('attr')
        def first(self):
            return f"first {self._attr}"

        @property
        @cached_property_based_on('attr')
        def second(self):
            return f"second {self._attr}"

    foo = Foo(1)
    # Cache both properties
    assert foo.first == "first 1"
    assert foo.second == "second 1"

    foo._attr = 2

    # Reading `first`` recomputes it and records the new attribute value.
    # Reading `second` must recompute it as well, since attr value has changed.
    assert foo.first == "first 2"
    assert foo.second == "second 2"


def test_cached_property_based_on_plain_attribute():
    class Foo:
        def __init__(self, attr):
            self.attr = attr

        @property
        @cached_property_based_on('attr')
        def prop(self):
            return self.attr * 10

    foo = Foo(1)
    assert foo.prop == 10
    foo.attr = 2
    assert foo.attr == 2
    assert foo.prop == 20


def test_cached_property_based_on_none_always_recomputes():
    """
    Regression test for https://github.com/sunpy/sunpy/issues/8780

    If the attribute that cache-invalidation is based on evaluates to
    `None` (e.g. because computing it failed), the property must always be
    recomputed, since `None == None` cannot be used to conclude that
    nothing has changed.
    """
    class Foo:
        def __init__(self, value):
            self._value = value

        @property
        def attr(self):
            # Always `None`, simulating an attribute whose value cannot be
            # reliably computed (e.g. `MetaDict.item_hash` returning `None`
            # because the metadata contains an unhashable value).
            return None

        @property
        @cached_property_based_on('attr')
        def prop(self):
            return self._value

    foo = Foo(1)
    assert foo.prop == 1

    # Even though `attr` is unchanged (`None`), the property must be
    # recomputed on every access, because a `None` attribute value means
    # "could not determine whether anything changed". This is asserted by
    # changing the underlying value and checking that `prop` picks it up,
    # which would not happen if a stale cached value were returned.
    foo._value = 2
    assert foo.prop == 2
    foo._value = 3
    assert foo.prop == 3
