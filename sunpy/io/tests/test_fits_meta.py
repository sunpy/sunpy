import pytest

import fits_meta


def test_SolarHeader_class_basic_creation():
    """
    """
    header = fits_meta.SolarnetHeader()
    assert isinstance(header, fits_meta.SolarnetHeader)