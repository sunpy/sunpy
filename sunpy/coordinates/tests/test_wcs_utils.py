
import numpy as np
import pytest

import astropy.units as u
from astropy.coordinates import BaseCoordinateFrame
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates.frames import (
    Heliocentric,
    HeliographicCarrington,
    HeliographicStonyhurst,
    Helioprojective,
)
from sunpy.util import SunpyUserWarning
from ..wcs_utils import _set_wcs_aux_obs_coord, solar_frame_to_wcs_mapping, solar_wcs_frame_mapping


@pytest.mark.parametrize('ctype, frame', [[['HPLN', 'HPLT'], Helioprojective],
                                          [['HPLT', 'HPLN'], Helioprojective],
                                          [['HGLN', 'HGLT'], HeliographicStonyhurst],
                                          [['CRLN', 'CRLT'], HeliographicCarrington],
                                          [['SOLX', 'SOLY'], Heliocentric]
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
    Make sure auxillary information round trips properly from coordinate frames
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
    assert result.rsun.value == header['rsun_ref']


def test_hpc_frame_to_wcs():
    frame = Helioprojective(observer="earth", obstime='2013-10-28')
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'HPLN-TAN'
    assert result_wcs.wcs.cunit[0] == 'arcsec'
    assert result_wcs.wcs.dateobs == '2013-10-28T00:00:00.000'

    new_frame = solar_wcs_frame_mapping(result_wcs)
    assert isinstance(new_frame.observer, HeliographicStonyhurst)
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
    frame = HeliographicStonyhurst(obstime='2013-10-28')
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'HGLN-TAN'
    assert result_wcs.wcs.cunit[0] == 'deg'
    assert result_wcs.wcs.dateobs == '2013-10-28T00:00:00.000'

    # Test a frame with no obstime
    frame = HeliographicStonyhurst()
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'HGLN-TAN'
    assert result_wcs.wcs.cunit[0] == 'deg'
    assert result_wcs.wcs.dateobs == ''


def test_hgc_frame_to_wcs():
    frame = HeliographicCarrington(obstime='2013-10-28')
    result_wcs = solar_frame_to_wcs_mapping(frame)

    assert isinstance(result_wcs, WCS)

    assert result_wcs.wcs.ctype[0] == 'CRLN-TAN'
    assert result_wcs.wcs.cunit[0] == 'deg'
    assert result_wcs.wcs.dateobs == '2013-10-28T00:00:00.000'

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


def test_attribute_warnings():
    # Check that warnings are raised if we try to convert a WCS with .rsun
    # or .heliographic_observer attributes
    wcs = WCS(naxis=2)
    wcs.rsun = None
    with pytest.warns(SunpyUserWarning,
                      match='Support for the .rsun attribute on a WCS is deprecated'):
        solar_wcs_frame_mapping(wcs)

    wcs = WCS(naxis=2)
    wcs.heliographic_observer = None
    with pytest.warns(SunpyUserWarning,
                      match='upport for the .heliographic_observer attribute on a WCS is deprecated'):
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
