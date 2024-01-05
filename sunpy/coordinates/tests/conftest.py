import pytest

from astropy.coordinates import solar_system_ephemeris


@pytest.fixture
def use_DE440s():
    # This fixture is for test functions that want to use the JPL DE440s ephemeris
    old_ephemeris = solar_system_ephemeris.get()
    solar_system_ephemeris.set('de440s')
    yield
    solar_system_ephemeris.set(old_ephemeris)
