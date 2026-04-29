import astropy.units as u
import pytest
import sunpy.net.attrs as a
from sunpy.util.exceptions import SunpyUserWarning

from sunpy.net.soar._attrs import SOOP, Distance, Product, Sensor, walker


def _apply(attr):
    """Apply a single attr via the walker and return the resulting params list."""
    params = []
    walker.apply(attr, params)
    return params


def test_time():
    attr = a.Time("2023-01-15 06:30:00", "2023-01-16 12:45:00")
    params = _apply(attr)
    assert len(params) == 1
    assert params[0] == "begin_time>='2023-01-15 06:30:00' AND begin_time<='2023-01-16 12:45:00'"


def test_level_from_int():
    """Integer level should be converted to 'L<n>' string."""
    params = _apply(a.Level(2))
    assert params == ["level='L2'"]


def test_level_from_string():
    """String levels should be uppercased."""
    params = _apply(a.Level("l1"))
    assert params == ["level='L1'"]


def test_level_low_latency():
    """Low-latency level strings should pass through uppercased."""
    params = _apply(a.Level("LL02"))
    assert params == ["level='LL02'"]


def test_level_invalid_warns():
    """An unrecognised level should emit a SunpyUserWarning."""
    with pytest.warns(SunpyUserWarning, match="level not in list of allowed levels"):
        _apply(a.Level("L9"))


def test_instrument():
    params = _apply(a.Instrument("EUI"))
    assert params == ["instrument='EUI'"]


def test_product():
    params = _apply(Product("eui-fsi174-image"))
    assert params == ["descriptor='eui-fsi174-image'"]


def test_provider():
    params = _apply(a.Provider("SOAR"))
    assert params == ["provider='SOAR'"]


def test_soop():
    params = _apply(SOOP("r_small_mres_mcad_ar_long_term"))
    assert params == ["soop_name='r_small_mres_mcad_ar_long_term'"]


def test_detector():
    params = _apply(a.Detector("HRI_EUV"))
    assert params == ["Detector='HRI_EUV'"]


def test_sensor():
    params = _apply(Sensor("ept"))
    assert params == ["Sensor='ept'"]


def test_wavelength():
    params = _apply(a.Wavelength(304 * u.AA, 304 * u.AA))
    assert len(params) == 1
    assert params[0] == "Wavemin='304.0' AND Wavemax='304.0'"


def test_wavelength_range():
    params = _apply(a.Wavelength(171 * u.AA, 304 * u.AA))
    assert params[0] == "Wavemin='171.0' AND Wavemax='304.0'"


def test_distance():
    params = _apply(Distance(0.3 * u.AU, 0.5 * u.AU))
    assert len(params) == 1
    assert params[0] == "DISTANCE(0.3,0.5)"


def test_distance_out_of_bounds_warns():
    """Distance values outside 0.28-1.0 AU should emit a warning."""
    with pytest.warns(SunpyUserWarning, match="Distance values must be within the range"):
        _apply(Distance(0.1 * u.AU, 0.5 * u.AU))


def test_distance_max_out_of_bounds_warns():
    """Only the max being out of bounds should also warn."""
    with pytest.warns(SunpyUserWarning, match="Distance values must be within the range"):
        _apply(Distance(0.3 * u.AU, 1.5 * u.AU))


def test_and_applier():
    """AttrAnd applier should recurse and collect all sub-attr params."""
    compound = a.Instrument("EUI") & a.Level(1) & Product("eui-fsi174-image")
    params = _apply(compound)
    assert "instrument='EUI'" in params
    assert "level='L1'" in params
    assert "descriptor='eui-fsi174-image'" in params
    assert len(params) == 3


def test_create_from_single_attr():
    """A single DataAttr should produce a list containing one sub-list."""
    result = walker.create(a.Instrument("EUI"))
    assert len(result) == 1
    assert result == [["instrument='EUI'"]]


def test_create_from_and():
    """An AND of attrs should produce a single sub-list with all params."""
    query = a.Instrument("EUI") & a.Level(2)
    result = walker.create(query)
    assert len(result) == 1
    assert "instrument='EUI'" in result[0]
    assert "level='L2'" in result[0]


def test_create_from_or():
    """An OR should produce multiple sub-lists, one per branch."""
    query = a.Instrument("EUI") | a.Instrument("STIX")
    result = walker.create(query)
    assert len(result) == 2
    assert result[0] == [["instrument='EUI'"]]
    assert result[1] == [["instrument='STIX'"]]


def test_create_from_compound_or_and():
    """OR of two AND branches should produce two sub-lists."""
    query = (a.Instrument("EUI") & a.Level(1)) | (a.Instrument("STIX") & a.Level(2))
    result = walker.create(query)
    assert len(result) == 2
    # Each branch is a list of one sub-list
    eui_params = result[0][0]
    stix_params = result[1][0]
    assert "instrument='EUI'" in eui_params
    assert "level='L1'" in eui_params
    assert "instrument='STIX'" in stix_params
    assert "level='L2'" in stix_params


def test_lowercase():
    """Product value should always be lowercased."""
    p = Product("EUI-FSI174-IMAGE")
    assert p.value == "eui-fsi174-image"


def test_already_lowercase():
    p = Product("eui-fsi174-image")
    assert p.value == "eui-fsi174-image"


def test_mixed_case():
    p = Product("Eui-Fsi174-Image")
    assert p.value == "eui-fsi174-image"


def test_converts_km_to_au():
    """Distance should convert km inputs to AU."""
    d = Distance(0.5 * u.AU.to(u.km) * u.km, 0.6 * u.AU.to(u.km) * u.km)
    assert u.isclose(d.min, 0.5 * u.AU)
    assert u.isclose(d.max, 0.6 * u.AU)


def test_au_passthrough():
    """AU inputs should be stored directly."""
    d = Distance(0.3 * u.AU, 0.4 * u.AU)
    assert u.isclose(d.min, 0.3 * u.AU)
    assert u.isclose(d.max, 0.4 * u.AU)


def test_non_scalar_raises():
    """Non-scalar quantities should raise ValueError."""
    with pytest.raises(ValueError, match="scalar"):
        Distance([0.3, 0.4] * u.AU, 0.5 * u.AU)


def test_non_scalar_max_raises():
    with pytest.raises(ValueError, match="scalar"):
        Distance(0.3 * u.AU, [0.4, 0.5] * u.AU)


def test_collides_with_same_class():
    d1 = Distance(0.3 * u.AU, 0.4 * u.AU)
    d2 = Distance(0.5 * u.AU, 0.6 * u.AU)
    assert d1.collides(d2) is True


def test_does_not_collide_with_other_type():
    d = Distance(0.3 * u.AU, 0.4 * u.AU)
    assert d.collides(a.Time("2023-01-01", "2023-01-02")) is False
