import warnings

import pytest

from sunpy.util.decorators import deprecated, get_removal_version
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
