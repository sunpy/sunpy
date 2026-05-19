import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from astropy.utils.introspection import minversion

from sunpy.coordinates.frames import HeliographicStonyhurst, Helioprojective


@pytest.fixture(params=[True, False], ids=["single_observer", "array_observer"])
def coords(request):
    single_observer = request.param
    times = Time("2025-01-01T00:00:00") + np.arange(11) * u.day

    if single_observer:
        observer = HeliographicStonyhurst(
            10 * u.deg,
            -5 * u.deg,
            0.8 * u.AU,
            obstime=times[0],
        )
        frame = Helioprojective(obstime=times[0], observer=observer)
    else:
        observer = HeliographicStonyhurst(
            np.linspace(10, 20, len(times)) * u.deg,
            np.linspace(-5, 5, len(times)) * u.deg,
            u.Quantity(np.full(len(times), 0.8), u.AU),
            obstime=times,
        )
        frame = Helioprojective(obstime=times, observer=observer)

    return SkyCoord(
        np.arange(len(times)) * u.arcsec,
        -1 * np.arange(len(times)) * u.arcsec,
        frame=frame,
    )


@pytest.mark.parametrize("fmt", ["fits", "ascii.ecsv", "parquet"])
def test_qtable_sunpy_coordinate_roundtrip(tmp_path, coords, fmt):
    # Older versions of astropy (< 7.0) do not support coordinate frames as mixins in tables,
    # specifically for array-valued attributes like the observer in this test.
    if not coords.observer.isscalar and not minversion("astropy", "7.0"):
        pytest.skip("Array-valued frame attributes require astropy >= 7.0 for serialization in some formats.")

    original = QTable()
    original["pos"] = coords

    output_file = tmp_path / f"sunpy_coords.{fmt}"
    original.write(output_file, format=fmt)

    roundtrip = QTable.read(output_file, format=fmt)
    result = roundtrip["pos"]

    assert isinstance(result.frame, Helioprojective)
    assert isinstance(result.observer, HeliographicStonyhurst)
    assert result.observer.isscalar is coords.observer.isscalar

    assert_quantity_allclose(original["pos"].Tx, result.Tx)
    assert_quantity_allclose(original["pos"].Ty, result.Ty)
    assert_quantity_allclose(original["pos"].observer.lon, result.observer.lon)
    assert_quantity_allclose(original["pos"].observer.lat, result.observer.lat)
    assert_quantity_allclose(original["pos"].observer.radius, result.observer.radius)
