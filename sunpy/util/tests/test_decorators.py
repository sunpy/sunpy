import warnings

import pytest

from sunpy.util.decorators import ACTIVE_CONTEXTS, deprecated, sunpycontextmanager
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
