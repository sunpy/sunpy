import warnings

import pytest

from sunpy.util.decorators import _active_contexts, deprecated, sunpycontextmanager
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
    assert _active_contexts.stack == []

    with ctx1():
        # Check that the context is active while inside
        assert _active_contexts.stack == [ctx1_name]

        with ctx2():
            # Check nesting of contexts
            assert _active_contexts.stack == [ctx1_name, ctx2_name]

            with ctx1():
                # Check a repeated context in the nesting
                assert _active_contexts.stack == [ctx1_name, ctx2_name, ctx1_name]

            # Check that only the last context is removed and not its duplicate
            assert _active_contexts.stack == [ctx1_name, ctx2_name]

        assert _active_contexts.stack == [ctx1_name]

    assert _active_contexts.stack == []
