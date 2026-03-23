import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from sunpy.coordinates.frames import HeliographicStonyhurst, Helioprojective


def _make_coords(single_observer):
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


@pytest.mark.parametrize(
    "single_observer", [True, False], ids=["single_observer", "array_observer"]
)
def test_qtable_sunpy_coordinate_fits_roundtrip(tmp_path, single_observer):
    original = QTable()
    original["pos"] = _make_coords(single_observer)

    output_file = tmp_path / "sunpy_coords.fits"
    original.write(output_file)

    roundtrip = QTable.read(output_file)
    result = roundtrip["pos"]

    assert isinstance(result.frame, Helioprojective)
    assert isinstance(result.observer, HeliographicStonyhurst)
    assert result.observer.isscalar is single_observer

    assert_quantity_allclose(original["pos"].Tx, result.Tx)
    assert_quantity_allclose(original["pos"].Ty, result.Ty)
    assert_quantity_allclose(original["pos"].observer.lon, result.observer.lon)
    assert_quantity_allclose(original["pos"].observer.lat, result.observer.lat)
    assert_quantity_allclose(original["pos"].observer.radius, result.observer.radius)
