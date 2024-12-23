
import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame, SkyCoord
from astropy.coordinates.earth import EarthLocation
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates.frames import (
    Heliocentric,
    HeliographicCarrington,
    HeliographicStonyhurst,
    Helioprojective,
    SunPyBaseCoordinateFrame,
)
from sunpy.coordinates.wcs_utils import (
    _set_wcs_aux_obs_coord,
    solar_frame_to_wcs_mapping,
    solar_wcs_frame_mapping,
)


@pytest.mark.parametrize(('ctype', 'frame'), [(['HPLN', 'HPLT'], Helioprojective),
                                              (['HPLT', 'HPLN'], Helioprojective),
                                              (['HGLN', 'HGLT'], HeliographicStonyhurst),
                                              (['CRLN', 'CRLT'], HeliographicCarrington),
                                              (['SOLX', 'SOLY'], Heliocentric)
                                              ])
def test_wcs_frame_mapping(ctype, frame):
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ctype
    result = solar_wcs_frame_mapping(wcs)
    assert isinstance(result, frame)


def test_wcs_frame_mapping_none():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['spam', 'eggs']

    result = solar_wcs_frame_mapping(wcs)

    assert result is None


def test_wcs_frame_mapping_dateavg():
    wcs = WCS(naxis=2)
    wcs.wcs.dateavg = "2020-01-01T00:00:00"
    wcs.wcs.ctype = ["HPLN", "HPLT"]

    result = solar_wcs_frame_mapping(wcs)

    assert result.obstime == Time("2020-01-01T00:00:00")

    wcs = WCS(naxis=2)
    wcs.wcs.dateavg = "2020-01-01T00:00:00"
    wcs.wcs.dateobs = "2020-01-01T01:00:00"
    wcs.wcs.ctype = ["HPLN", "HPLT"]

    result = solar_wcs_frame_mapping(wcs)

    assert result.obstime == Time("2020-01-01T00:00:00")


def test_wcs_frame_mapping_observer_hgc_self():
    # Test whether a WCS with HGC coordinates for the observer location uses observer="self"
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['SOLX', 'SOLY']
    wcs.wcs.dateobs = '2001-01-01'
    wcs.wcs.aux.crln_obs = 10
    wcs.wcs.aux.hglt_obs = 20
    wcs.wcs.aux.dsun_obs = 1.5e11

    # This frame will have the observer location in HGS
    result = solar_wcs_frame_mapping(wcs)

    # We perform the desired transformation using observer="self"
    hgc_obs = HeliographicCarrington(wcs.wcs.aux.crln_obs * u.deg,
                                     wcs.wcs.aux.hglt_obs * u.deg,
                                     wcs.wcs.aux.dsun_obs * u.m,
                                     obstime=wcs.wcs.dateobs, observer="self")
    hgs_obs = hgc_obs.transform_to(HeliographicStonyhurst(obstime=hgc_obs.obstime))

    assert_quantity_allclose(result.observer.lon, hgs_obs.lon)
    assert_quantity_allclose(result.observer.lat, hgs_obs.lat)
    assert_quantity_allclose(result.observer.radius, hgs_obs.radius)


def test_wcs_aux():
    """
    Make sure auxiliary information round trips properly from coordinate frames
    to WCS and back.
    """
    data = np.ones([6, 6], dtype=np.float64)
    header = {'CRVAL1': 0,
              'CRVAL2': 0,
              'CRPIX1': 5,
              'CRPIX2': 5,
              'CDELT1': 10,
              'CDELT2': 10,
              'CUNIT1': 'arcsec',
              'CUNIT2': 'arcsec',
              'PC1_1': 0,
              'PC1_2': -1,
              'PC2_1': 1,
              'PC2_2': 0,
              'NAXIS1': 6,
              'NAXIS2': 6,
              'CTYPE1': 'HPLN-TAN',
              'CTYPE2': 'HPLT-TAN',
              'date-obs': '1970-01-01T00:00:00',
              'obsrvtry': 'Foo',
              'detector': 'bar',
              'wavelnth': 10,
              'waveunit': 'm',
              'hglt_obs': 0,
              'hgln_obs': 0,
              'dsun_obs': 10,
              'rsun_ref': 690000000}
    generic_map = sunpy.map.Map((data, header))

    wcs = generic_map.wcs
    assert wcs.wcs.aux.hglt_obs == 0
    assert wcs.wcs.aux.hgln_obs == 0
    assert wcs.wcs.aux.dsun_obs == 10
    assert wcs.wcs.aux.rsun_ref == header['rsun_ref']

    result = solar_wcs_frame_mapping(wcs)

    assert isinstance(result, Helioprojective)
    assert result.observer.lat.value == 0
    assert result.observer.lon.value == 0
    assert result.observer.radius.value == 10
    assert result.observer.rsun.to_value(u.m) == header['rsun_ref']
    assert result.rsun.to_value(u.m) == header['rsun_ref']


def test_hpc_frame_to_wcs():
    frame = Helioprojective(observer="earth", obstime='2013-10-28', rsun=690*u.Mm)
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'HPLN-TAN'
    assert result_wcs.wcs.cunit[0] == 'arcsec'
    assert result_wcs.wcs.dateobs == '2013-10-28T00:00:00.000'
    assert result_wcs.wcs.aux.rsun_ref == frame.rsun.to_value(u.m)

    new_frame = solar_wcs_frame_mapping(result_wcs)
    assert new_frame.rsun == frame.rsun

    # Test a frame with no obstime and no observer
    frame = Helioprojective()
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'HPLN-TAN'
    assert result_wcs.wcs.cunit[0] == 'arcsec'
    assert result_wcs.wcs.dateobs == ''

    new_frame = solar_wcs_frame_mapping(result_wcs)
    assert new_frame.observer is None
    assert new_frame.rsun == frame.rsun


def test_hgs_frame_to_wcs():
    frame = HeliographicStonyhurst(obstime='2013-10-28', rsun=690*u.Mm)
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'HGLN-TAN'
    assert result_wcs.wcs.cunit[0] == 'deg'
    assert result_wcs.wcs.dateobs == '2013-10-28T00:00:00.000'
    assert result_wcs.wcs.aux.rsun_ref == frame.rsun.to_value(u.m)

    new_frame = solar_wcs_frame_mapping(result_wcs)
    assert new_frame.rsun == frame.rsun

    # Test a frame with no obstime
    frame = HeliographicStonyhurst()
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'HGLN-TAN'
    assert result_wcs.wcs.cunit[0] == 'deg'
    assert result_wcs.wcs.dateobs == ''


def test_hgc_frame_to_wcs():
    frame = HeliographicCarrington(obstime='2013-10-28', rsun=690*u.Mm)
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'CRLN-TAN'
    assert result_wcs.wcs.cunit[0] == 'deg'
    assert result_wcs.wcs.dateobs == '2013-10-28T00:00:00.000'
    assert result_wcs.wcs.aux.rsun_ref == frame.rsun.to_value(u.m)

    new_frame = solar_wcs_frame_mapping(result_wcs)
    assert new_frame.rsun == frame.rsun

    # Test a frame with no obstime
    frame = HeliographicCarrington()
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'CRLN-TAN'
    assert result_wcs.wcs.cunit[0] == 'deg'
    assert result_wcs.wcs.dateobs == ''


def test_hcc_frame_to_wcs():
    frame = Heliocentric(observer="earth", obstime='2013-10-28')
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'SOLX'
    assert result_wcs.wcs.dateobs == '2013-10-28T00:00:00.000'

    new_frame = solar_wcs_frame_mapping(result_wcs)
    assert isinstance(new_frame.observer, HeliographicStonyhurst)

    # Test a frame with no obstime and no observer
    frame = Heliocentric()
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'SOLX'
    assert result_wcs.wcs.dateobs == ''

    new_frame = solar_wcs_frame_mapping(result_wcs)
    assert new_frame.observer is None


def test_non_sunpy_frame_to_wcs():
    # For a non-SunPy frame, our mapping should return None
    frame = BaseCoordinateFrame()
    assert solar_frame_to_wcs_mapping(frame) is None


def test_attribute_errors():
    # Check that errors are raised if we try to convert a WCS with .rsun
    # or .heliographic_observer attributes
    wcs = WCS(naxis=2)
    wcs.rsun = None
    with pytest.raises(ValueError, match='The .rsun attribute'):
        solar_wcs_frame_mapping(wcs)

    wcs = WCS(naxis=2)
    wcs.heliographic_observer = None
    with pytest.raises(ValueError, match='The .heliographic_observer attribute'):
        solar_wcs_frame_mapping(wcs)

    wcs = WCS(naxis=2)
    wcs.heliographic_observer = None
    wcs.rsun = None
    with pytest.raises(ValueError, match='The .rsun and .heliographic_observer attribute'):
        solar_wcs_frame_mapping(wcs)


def test_set_wcs_aux():
    wcs = WCS(naxis=2)
    observer = Helioprojective(observer="earth", obstime='2013-10-28').observer
    _set_wcs_aux_obs_coord(wcs, observer)
    assert wcs.wcs.aux.hgln_obs == 0
    assert u.allclose(wcs.wcs.aux.hglt_obs, 4.7711570596394015)
    assert u.allclose(wcs.wcs.aux.dsun_obs, 148644585949.4918)
    assert wcs.wcs.aux.crln_obs is None

    wcs = WCS(naxis=2)
    observer = observer.transform_to(HeliographicCarrington(observer=observer))
    _set_wcs_aux_obs_coord(wcs, observer)
    assert wcs.wcs.aux.hgln_obs is None
    assert u.allclose(wcs.wcs.aux.hglt_obs, 4.7711570596394015)
    assert u.allclose(wcs.wcs.aux.dsun_obs, 148644585949.4918)
    assert u.allclose(wcs.wcs.aux.crln_obs, 326.05139910339886)

    observer = observer.transform_to(Heliocentric(observer=observer))
    with pytest.raises(ValueError, match='obs_coord must be in a Stonyhurst or Carrington frame'):
        _set_wcs_aux_obs_coord(wcs, observer)


def test_self_observer():
    frame = HeliographicCarrington(0*u.deg, 0*u.deg, 1*u.au,
                                   observer="self", obstime='2013-10-28')
    wcs = solar_frame_to_wcs_mapping(frame)
    assert wcs.wcs.aux.hgln_obs is None
    assert u.allclose(wcs.wcs.aux.hglt_obs, frame.lon.to_value(u.deg))
    assert u.allclose(wcs.wcs.aux.crln_obs, frame.lon.to_value(u.deg))
    assert u.allclose(wcs.wcs.aux.dsun_obs, frame.radius.to_value(u.m))


@pytest.fixture
def dkist_location():
    return EarthLocation(*(-5466045.25695494, -2404388.73741278, 2242133.88769004) * u.m)


@pytest.mark.remote_data
def test_obsgeo_frame_mapping_cartesian(dkist_location, caplog):

    obstime = Time("2021-05-21T03:00:00")
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['HPLT', 'HPLN']
    wcs.wcs.obsgeo = list(dkist_location.to_value(u.m).tolist()) + [0, 0, 0]
    wcs.wcs.dateobs = obstime.isot

    frame = solar_wcs_frame_mapping(wcs)

    assert frame.observer is not None

    assert frame.observer == SkyCoord(dkist_location.get_itrs(obstime)).transform_to('heliographic_stonyhurst').frame

    assert not caplog.records


@pytest.mark.remote_data
def test_frame_mapping_obsgeo_spherical(dkist_location, caplog):

    obstime = Time("2021-05-21T03:00:00")
    location = dkist_location.get_itrs(obstime)
    loc_sph = location.spherical
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ['HPLT', 'HPLN']
    wcs.wcs.obsgeo = [0, 0, 0] + [loc_sph.lon.to_value(u.deg), loc_sph.lat.to_value(u.deg), loc_sph.distance.to_value(u.m)]
    wcs.wcs.dateobs = obstime.isot

    frame = solar_wcs_frame_mapping(wcs)

    assert frame.observer is not None

    assert frame.observer == SkyCoord(location).transform_to('heliographic_stonyhurst').frame

    assert not caplog.records


def test_observer_hgln_crln_priority():
    """
    When extracting the observer information from a FITS header, ensure
    Stonyhurst (HG) coordinates are preferred over Carrington (CR) if present.
    Note that `Map` creates a custom FITS header with a sanitized observer
    location, so it is separately tested in the map module.
    """
    header = {'CRVAL1': 0,
              'CRVAL2': 0,
              'CRPIX1': 5,
              'CRPIX2': 5,
              'CDELT1': 10,
              'CDELT2': 10,
              'CUNIT1': 'arcsec',
              'CUNIT2': 'arcsec',
              'PC1_1': 0,
              'PC1_2': -1,
              'PC2_1': 1,
              'PC2_2': 0,
              'NAXIS1': 6,
              'NAXIS2': 6,
              'CTYPE1': 'HPLN-TAN',
              'CTYPE2': 'HPLT-TAN',
              'date-obs': '1970-01-01T00:00:00',
              'mjd-obs': 40587,
              'hglt_obs': 0,
              'hgln_obs': 0,
              'crlt_obs': 2,
              'crln_obs': 2,
              'dsun_obs': 10,
              'rsun_ref': 690000000}
    wcs = WCS(header)
    c = wcs.pixel_to_world(0, 0)
    assert c.observer.lon == 0 * u.deg
    # Note: don't test whether crlt or hglt is used---according to
    # _set_wcs_aux_obs_coord, those are expected to always be the same and so
    # the same one is always used


def test_sunpybaseframe_external():
    class MyFrame(SunPyBaseCoordinateFrame):
        pass

    out = solar_frame_to_wcs_mapping(MyFrame())
    assert out is None
