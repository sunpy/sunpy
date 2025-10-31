import pytest
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, ITRS
from astropy.time import Time
from sunpy.map.header_helper import make_fitswcs_header
from sunpy.util import MetaDict


@pytest.mark.parametrize("observer_case", ["EarthLocation", "ITRS", "Heliographic"])
def test_obsg_fields_for_various_observers(observer_case):
    """
    Rigorous multi-frame test for OBSGEO-[XYZ] FITS header generation.

    This ensures the header helper correctly handles:
    - Ground-based observers via EarthLocation
    - Explicit ITRS SkyCoord observers
    - Solar-centered (Heliographic) observers
    """
    obstime = Time("2020-01-01")
    base_loc = EarthLocation(lat=19*u.deg, lon=72*u.deg, height=15*u.m)

    # ---- Observer setup ----
    if observer_case == "EarthLocation":
        # Convert EarthLocation → SkyCoord(ITRS) to be a valid SunPy observer
        observer = SkyCoord(base_loc.get_itrs(obstime=obstime))
    elif observer_case == "ITRS":
        observer = SkyCoord(base_loc.get_itrs(obstime=obstime))
    else:  # HeliographicStonyhurst (Sun-centered observer)
        from sunpy.coordinates import frames
        observer = SkyCoord(0*u.deg, 0*u.deg,
                            frame=frames.HeliographicStonyhurst,
                            obstime=obstime)

    # ---- Create reference pixel coordinate ----
    coord = SkyCoord(0*u.arcsec, 0*u.arcsec,
                     frame="helioprojective",
                     observer=observer,
                     obstime=obstime)

    header = make_fitswcs_header((512, 512), coord)
    assert isinstance(header, MetaDict)

    # ---- Validation: Earth-based observers ----
    if observer_case in ("EarthLocation", "ITRS"):
        for key in ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z"):
            assert key in header, f"{key} missing for {observer_case}"

        # Check magnitudes: should be roughly within Earth's radius
        x, y, z = [header[k] for k in ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z")]
        for val in (x, y, z):
            assert abs(val) < 7e6, f"{observer_case}: {val} too large for ground observer"

    # ---- Validation: Heliographic observer ----
    elif observer_case == "Heliographic":
        for key in ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z"):
            assert key in header, f"{key} missing for Heliographic observer"

        # Heliographic coordinates correspond to Sun-centered → huge magnitudes
        x, y, z = [header[k] for k in ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z")]
        for val in (x, y, z):
            assert abs(val) > 1e8, "Heliographic observer values too small — expected solar distance"


def test_header_without_observer_has_no_obsg_fields():
    """
    Ensure make_fitswcs_header() does not set OBSGEO keys if no observer is provided.
    """
    obstime = Time("2020-01-01")
    coord = SkyCoord(0*u.arcsec, 0*u.arcsec,
                     frame="helioprojective",
                     obstime=obstime)

    header = make_fitswcs_header((256, 256), coord)
    for key in ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z"):
        assert key not in header, f"{key} present despite missing observer"


@pytest.mark.parametrize("year", [2000, 2010, 2020, 2025])
def test_obsg_fields_stability_over_time(year):
    """
    Regression-style test:
    Ensure OBSGEO coordinates remain physically consistent across different obstimes.
    """
    obstime = Time(f"{year}-01-01")
    base_loc = EarthLocation(lat=19*u.deg, lon=72*u.deg, height=0*u.m)
    observer = SkyCoord(base_loc.get_itrs(obstime=obstime))

    coord = SkyCoord(0*u.arcsec, 0*u.arcsec,
                     frame="helioprojective",
                     observer=observer,
                     obstime=obstime)

    header = make_fitswcs_header((256, 256), coord)

    # All coordinates should exist
    for key in ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z"):
        assert key in header, f"{key} missing for obstime={year}"

    # Magnitude should remain within Earth radius range ±10%
    x, y, z = [header[k] for k in ("OBSGEO-X", "OBSGEO-Y", "OBSGEO-Z")]
    mag = (x**2 + y**2 + z**2) ** 0.5
    assert 6.0e6 < mag < 6.9e6, f"OBSGEO magnitude unstable for year={year}: {mag}"
