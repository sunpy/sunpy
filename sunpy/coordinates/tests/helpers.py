import astropy.units as u
from astropy.coordinates import Longitude
from astropy.tests.helper import assert_quantity_allclose


def assert_longitude_allclose(actual, desired, atol=None):
    """
    This works like :func:`~astropy.tests.helper.assert_quantity_allclose`, except it is intended
    for longitude angles (i.e., angles that wrap every 360 degrees).
    """
    difference = Longitude(actual - desired, wrap_angle=180*u.deg)
    assert_quantity_allclose(difference, 0*u.deg, atol=atol)
