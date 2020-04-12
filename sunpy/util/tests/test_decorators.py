from sunpy.util.decorators import get_removal_version


def test_removal_version_since_lts():
    major, minor = get_removal_version('2.0')
    assert major == 2
    assert minor == 1


def test_removal_version_not_since_lts():
    major, minor = get_removal_version('2.1')
    assert major == 3
    assert minor == 1
