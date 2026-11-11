import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import CartesianRepresentation, SkyCoord

import sunpy.net.jsoc as jsoc
import sunpy.net.jsoc.attrs as attrs
from sunpy.coordinates import HeliographicStonyhurst, Helioprojective, frames, get_earth
from sunpy.net import _attrs as core_attrs
from sunpy.net import attrs as a
from sunpy.net.attr import AttrAnd, AttrOr
from sunpy.util.exceptions import SunpyUserWarning


@pytest.mark.parametrize((("attr1", "attr2")),
                         [(attrs.Series('foo'), attrs.Series('boo')),
                          (attrs.Protocol('a1'), attrs.Protocol('a2')),
                          (attrs.Notify('email@somemail.com'),
                           attrs.Notify('someemail@somemail.com'))])
def test_and(attr1, attr2):
    with pytest.raises(TypeError, match=r"unsupported operand type\(s\) for \&"):
        attr1 & attr2


def test_basicquery():
    a1 = attrs.Series('foo')
    t1 = core_attrs.Time('2012/01/01', '2013/1/2')
    ans1 = jsoc.jsoc.and_(a1, t1)
    assert isinstance(ans1, AttrAnd)
    assert len(ans1.attrs) == 2


def test_mediumquery():
    a1 = attrs.Series('foo1')
    a2 = attrs.Series('foo2')
    t1 = core_attrs.Time('2012/01/01', '2013/1/2')
    ans1 = jsoc.jsoc.and_(a1 | a2, t1)
    assert isinstance(ans1, AttrOr)
    assert isinstance(ans1.attrs[0], AttrAnd)
    assert isinstance(ans1.attrs[1], AttrAnd)


def test_complexquery():
    a1 = attrs.Series('foo1')
    a2 = attrs.Series('foo2')
    t1 = core_attrs.Time('2012/01/01', '2013/1/2')
    t2 = core_attrs.Time('2012/01/01', '2013/1/3')
    ans1 = jsoc.jsoc.and_(a1 | a2, t1 | t2)
    assert isinstance(ans1.attrs[0], AttrOr)
    assert isinstance(ans1.attrs[0].attrs[0], AttrAnd)
    assert isinstance(ans1.attrs[0].attrs[1], AttrAnd)


def test_wavelength_error():
    with pytest.raises(TypeError):
        attrs.Wavelength('wobble')
    with pytest.raises(TypeError):
        attrs.Wavelength(3.24)
    with pytest.raises(TypeError):
        attrs.Wavelength((3, 3))


def test_wave_self():
    w1 = attrs.Wavelength(193*u.AA)
    assert jsoc.jsoc.and_(w1 | w1) is w1


def test_duplicate():
    w1 = attrs.Wavelength(193*u.AA)
    w2 = attrs.Wavelength(193*u.AA)
    assert jsoc.jsoc.and_(w1 | w2).min is w1.min


def test_random():
    w1 = attrs.Wavelength(193*u.AA)
    w2 = attrs.Series('spam')
    assert jsoc.jsoc.and_(w1 | w2) == AttrOr([w1, w2])


def test_empty_notify():
    with pytest.raises(ValueError, match="Notify attribute must contain an email address"):
        attrs.Notify(None)


def test_not_email_notify():
    with pytest.raises(ValueError, match="Notify attribute must contain an '@' symbol to be a valid email address"):
        attrs.Notify("someemailthatisntone")


def test_cutout_not_helioprojective():
    hpc = SkyCoord(500*u.arcsec, -200*u.arcsec,
                   obstime='2025-09-16', observer="earth", frame=frames.Helioprojective)
    # No error because helioprojective
    _  = attrs.Cutout(hpc, width=900*u.arcsec, height=900*u.arcsec)

    hpr = SkyCoord(500*u.arcsec, -200*u.arcsec,
                   obstime='2025-09-16', observer="earth", frame=frames.HelioprojectiveRadial)
    # Error because not helioprojective
    with pytest.raises(ValueError, match="`bottom_left` must be in the `Helioprojective` frame, but is instead in the `HelioprojectiveRadial` frame"):
        _ = attrs.Cutout(hpr, width=900*u.arcsec, height=900*u.arcsec)


def test_cutout_not_on_disk_when_tracking():
    bottom_left = SkyCoord(500*u.arcsec, -200*u.arcsec,
                           obstime='2025-09-16', observer="earth", frame=frames.Helioprojective)

    # No error because tracking is disabled
    cutout  = attrs.Cutout(bottom_left, width=900*u.arcsec, height=900*u.arcsec, tracking=False)
    assert cutout.value["x"] == 950
    assert cutout.value["y"] == 250

    # Error because tracking is enabled
    with pytest.raises(ValueError, match="Tracking is enabled, but the center of the cutout .* is not on the solar disk"):
        _  = attrs.Cutout(bottom_left, width=900*u.arcsec, height=900*u.arcsec, tracking=True)
OBSTIME = '2012-09-24T14:56:03'


@pytest.fixture
def earth_frame():
    return Helioprojective(obstime=OBSTIME, observer='earth')


@pytest.fixture
def offset_observer():
    """An observer offset from Earth center by a geosynchronous radius."""
    cart = get_earth(OBSTIME).cartesian
    offset = CartesianRepresentation(0 * u.km, 42164 * u.km, 0 * u.km)
    return SkyCoord(cart + offset, frame=HeliographicStonyhurst(obstime=OBSTIME))


def _cutout(frame, bl_xy=(-500, -275), tr_xy=(150, 375), **kwargs):
    bl = SkyCoord(bl_xy[0] * u.arcsec, bl_xy[1] * u.arcsec, frame=frame)
    tr = SkyCoord(tr_xy[0] * u.arcsec, tr_xy[1] * u.arcsec, frame=frame)
    return a.jsoc.Cutout(bl, tr, **kwargs)


def test_earth_observer_is_passthrough(earth_frame):
    """An Earth-observer coordinate needs no transformation and must not warn."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error", SunpyUserWarning)
        value = _cutout(earth_frame).value
    assert value['x'] == pytest.approx(-175.0)
    assert value['y'] == pytest.approx(50.0)
    assert value['width'] == pytest.approx(650.0)
    assert value['height'] == pytest.approx(650.0)
    assert value['locunits'] == 'arcsec'


def test_non_earth_observer_is_transformed(earth_frame, offset_observer):
    """A non-Earth observer changes the request and warns about the assumptions."""
    frame = Helioprojective(obstime=OBSTIME, observer=offset_observer)
    with pytest.warns(SunpyUserWarning, match="transformed to a Helioprojective frame"):
        shifted = _cutout(frame).value
    reference = _cutout(earth_frame).value
    assert shifted['x'] != reference['x']
    # The observer offset is small, so the shift should be sub-arcsecond
    assert abs(shifted['x'] - reference['x']) < 1.0
    assert abs(shifted['y'] - reference['y']) < 1.0


def test_transformed_box_contains_request(offset_observer):
    """The returned box is the bounding box of the transformed corners."""
    frame = Helioprojective(obstime=OBSTIME, observer=offset_observer)
    with pytest.warns(SunpyUserWarning):
        value = _cutout(frame).value
    assert value['width'] >= 650.0
    assert value['height'] >= 650.0


def test_off_limb_coordinates_use_spherical_screen(offset_observer):
    """Off-limb corners must not produce NaNs in the request."""
    frame = Helioprojective(obstime=OBSTIME, observer=offset_observer)
    with pytest.warns(SunpyUserWarning):
        value = _cutout(frame, bl_xy=(1050, -50), tr_xy=(1150, 50)).value
    assert np.isfinite(value['x'])
    assert np.isfinite(value['y'])
    assert np.isfinite(value['width'])
    assert np.isfinite(value['height'])


def test_distant_observer_changes_request_substantially():
    """A genuinely different vantage point produces a very different request."""
    observer = SkyCoord(35 * u.deg, 5 * u.deg, 0.55 * u.AU,
                        frame=HeliographicStonyhurst(obstime=OBSTIME))
    frame = Helioprojective(obstime=OBSTIME, observer=observer)
    with pytest.warns(SunpyUserWarning):
        value = _cutout(frame, bl_xy=(-200, -200), tr_xy=(200, 200)).value
    assert abs(value['x']) > 100


def test_missing_observer_raises():
    frame = Helioprojective(obstime=OBSTIME, observer=None)
    with pytest.raises(ValueError, match="must have `observer` set"):
        _cutout(frame)


def test_missing_obstime_raises():
    frame = Helioprojective(obstime=None, observer='earth')
    with pytest.raises(ValueError, match="must have `obstime` set"):
        _cutout(frame)


def test_non_helioprojective_raises():
    frame = HeliographicStonyhurst(obstime=OBSTIME)
    bl = SkyCoord(10 * u.deg, 10 * u.deg, frame=frame)
    tr = SkyCoord(20 * u.deg, 20 * u.deg, frame=frame)
    with pytest.raises(ValueError, match="must be in the `Helioprojective` frame"):
        a.jsoc.Cutout(bl, tr)


def test_tracking_check_uses_transformed_center(offset_observer):
    """Tracking validation happens after the transformation, not before."""
    frame = Helioprojective(obstime=OBSTIME, observer=offset_observer)
    with pytest.warns(SunpyUserWarning):
        value = _cutout(frame, tracking=True).value
    assert value['t'] == 0  # tracking enabled -> JSOC 't' is the negation


def test_tracking_off_disk_still_raises(offset_observer):
    frame = Helioprojective(obstime=OBSTIME, observer=offset_observer)
    with pytest.raises(ValueError, match="not on the solar disk"):
        with pytest.warns(SunpyUserWarning):
            _cutout(frame, bl_xy=(1050, -50), tr_xy=(1150, 50), tracking=True)
