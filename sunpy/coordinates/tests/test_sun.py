import re
import warnings

import numpy as np
import pytest
from erfa import ErfaWarning
from numpy.testing import assert_allclose

import astropy.units as u
from astropy.coordinates import Angle, EarthLocation
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time

from sunpy.coordinates import sun
from sunpy.sun.constants import radius
from .helpers import assert_longitude_allclose

# Ignore warnings that astropy throws when trying and failing to update ephemeris
pytestmark = pytest.mark.filterwarnings("ignore:failed to download")


@pytest.fixture
def t1():
    return Time('1992-10-13', scale='tdb')


@pytest.fixture
def t2():
    return Time('2013-06-17', scale='tt')


def test_true_longitude(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C16)
    assert_quantity_allclose(sun.true_longitude(t1), Angle('199d54m26.17s'), atol=0.005*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C12)
    assert_quantity_allclose(sun.true_longitude(t2), Angle('85d58m57.46s'), atol=0.005*u.arcsec)


def test_apparent_longitude(t1, t2):
    aberration_1au = Angle('20.496s')
    # Validate against a published value from the Astronomical Almanac (1992, C2+C16+B30)
    assert_quantity_allclose(sun.apparent_longitude(t1),
                             Angle('199d54m26.17s') + Angle('15.908s') - aberration_1au / 0.9976085,
                             atol=0.005*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C2+C12+B61)
    assert_quantity_allclose(sun.apparent_longitude(t2),
                             Angle('85d58m57.46') + Angle('12.0489s') - aberration_1au / 1.0159149,
                             atol=0.005*u.arcsec)


def test_true_latitude(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C16)
    # Disagreement with published precision may be due to changes in definitions after 1992
    assert_quantity_allclose(sun.true_latitude(t1), Angle('0.72s'), atol=0.05*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C12)
    assert_quantity_allclose(sun.true_latitude(t2), Angle('-0.42s'), atol=0.005*u.arcsec)


def test_apparent_latitude(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C2+C16)
    # Disagreement with published precision may be due to changes in definitions after 1992
    assert_quantity_allclose(sun.apparent_latitude(t1), Angle('0.72s'), atol=0.05*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C2+C12)
    assert_quantity_allclose(sun.apparent_latitude(t2), Angle('-0.42s'), atol=0.005*u.arcsec)


def test_angular_radius(t2):
    # Validate against a published value from the Astronomical Almanac (2013, C13)
    # The Astromomical Almanac uses a slightly different radius for the Sun (6.96e5 km)
    # The Astronomical Almanac also uses a small-angle approximation
    assert_quantity_allclose(sun.angular_radius(t2),
                             Angle('0d15m44.61s') / (6.96e5*u.km) * radius,  # scale to IAU radius
                             atol=0.005*u.arcsec)

    # Regression-only test
    assert_quantity_allclose(sun.angular_radius("2012/11/11"), 968.875*u.arcsec, atol=1e-3*u.arcsec)
    with pytest.warns(ErfaWarning):
        assert_quantity_allclose(sun.angular_radius("2043/03/01"),
                                 968.330*u.arcsec, atol=1e-3*u.arcsec)
    assert_quantity_allclose(sun.angular_radius("2001/07/21"), 944.042*u.arcsec, atol=1e-3*u.arcsec)


def test_mean_obliquity_of_ecliptic(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, B30)
    # Note that the publication date pre-dates the IAU 2006 definition of obliquity
    assert_quantity_allclose(sun.mean_obliquity_of_ecliptic(t1),
                             Angle('23d26m24.519s') + 0.308*u.arcsec, atol=0.05*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, B61)
    assert_quantity_allclose(sun.mean_obliquity_of_ecliptic(t2),
                             Angle('23d26m08.0875s') + 7.0153*u.arcsec, atol=0.00005*u.arcsec)


def test_true_rightascension():
    # Regression-only test
    assert_quantity_allclose(sun.true_rightascension("2012/11/11"), 226.548*u.deg, atol=1e-3*u.deg)
    with pytest.warns(ErfaWarning):
        assert_quantity_allclose(sun.true_rightascension(
            "2142/02/03"), 316.466*u.deg, atol=1e-3*u.deg)
    assert_quantity_allclose(sun.true_rightascension("2013/12/11"), 258.150*u.deg, atol=1e-3*u.deg)


def test_true_rightascension_J2000(t1, t2):
    # Validate against JPL HORIZONS output
    # https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE=OBSERVER&COMMAND=10&CENTER=500&QUANTITIES=1&ANG_FORMAT=HMS&EXTRA_PREC=YES&TIME_TYPE=TT&TLIST=2448908.5%0A2456460.5
    assert_quantity_allclose(sun.true_rightascension(t1, equinox_of_date=False),
                             Angle('13h13m53.65s'), atol=Angle('0h0m0.005s'))
    assert_quantity_allclose(sun.true_rightascension(t2, equinox_of_date=False),
                             Angle('05h41m40.32s'), atol=Angle('0h0m0.005s'))


def test_true_declination():
    # Regression-only test
    assert_quantity_allclose(sun.true_declination("2012/11/11"), -17.470*u.deg, atol=1e-3*u.deg)
    with pytest.warns(ErfaWarning):
        assert_quantity_allclose(sun.true_declination("2245/12/01"), -21.717*u.deg, atol=1e-3*u.deg)
    assert_quantity_allclose(sun.true_declination("2014/05/27"), 21.245*u.deg, atol=1e-3*u.deg)


def test_true_declination_J2000(t1, t2):
    # Validate against JPL HORIZONS output
    # https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE=OBSERVER&COMMAND=10&CENTER=500&QUANTITIES=1&ANG_FORMAT=HMS&EXTRA_PREC=YES&TIME_TYPE=TT&TLIST=2448908.5%0A2456460.5
    assert_quantity_allclose(sun.true_declination(t1, equinox_of_date=False),
                             Angle('-7d49m20.84s'), atol=0.005*u.arcsec)
    assert_quantity_allclose(sun.true_declination(t2, equinox_of_date=False),
                             Angle('23d22m13.97s'), atol=0.005*u.arcsec)


def test_true_obliquity_of_ecliptic(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, B30)
    # Note that the publication date pre-dates the IAU 2006 definition of obliquity
    assert_quantity_allclose(sun.true_obliquity_of_ecliptic(t1), Angle('23d26m24.519s'),
                             atol=0.05*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, B61)
    assert_quantity_allclose(sun.true_obliquity_of_ecliptic(t2), Angle('23d26m08.0875'),
                             atol=0.00005*u.arcsec)


def test_apparent_rightascension(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C16)
    assert_quantity_allclose(sun.apparent_rightascension(t1), Angle('13h13m30.75s'),
                             atol=Angle('0h0m0.005s'))

    # Validate against a published value from the Astronomical Almanac (2013, C12)
    assert_quantity_allclose(sun.apparent_rightascension(t2), Angle('5h42m28.88s'),
                             atol=Angle('0h0m0.005s'))


def test_apparent_rightascension_J2000(t1):
    # Regression-only test
    assert_quantity_allclose(sun.apparent_rightascension(t1, equinox_of_date=False),
                             Angle('13h13m52.37s'), atol=Angle('0h0m0.005s'))


def test_apparent_declination(t1, t2):
    # Validate against a published value from the Astronomical Almanac (1992, C16)
    assert_quantity_allclose(sun.apparent_declination(t1), Angle('-7d47m01.7s'), atol=0.05*u.arcsec)

    # Validate against a published value from the Astronomical Almanac (2013, C12)
    assert_quantity_allclose(sun.apparent_declination(t2), Angle('23d22m27.8s'), atol=0.05*u.arcsec)


def test_apparent_declination_J2000(t1):
    # Regression-only test
    assert_quantity_allclose(sun.apparent_declination(t1, equinox_of_date=False),
                             Angle('-7d49m13.1s'), atol=0.05*u.arcsec)


def test_sky_position(t1, t2):
    pos1 = sun.sky_position(t1)
    ra1 = sun.apparent_rightascension(t1)
    dec1 = sun.apparent_declination(t1)
    assert_quantity_allclose(pos1, (ra1, dec1))

    pos2 = sun.sky_position(t2, equinox_of_date=False)
    ra2 = sun.apparent_rightascension(t2, equinox_of_date=False)
    dec2 = sun.apparent_declination(t2, equinox_of_date=False)
    assert_quantity_allclose(pos2, (ra2, dec2))


@pytest.mark.filterwarnings('ignore:Tried to get polar motions for times after IERS data is valid')
@pytest.mark.filterwarnings(r'ignore:\(some\) times are outside of range covered by IERS table')
def test_print_params():
    # Test only for any issues with printing; accuracy is covered by other tests
    sun.print_params()


def test_B0():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(sun.B0('1992-Oct-13'), 5.99*u.deg, atol=5e-3*u.deg)


def test_B0_astronomical_almanac():
    # Validate against published values from the Astronomical Almanac (2013)
    sun_B0 = sun.B0(Time(['2013-04-01', '2013-12-01'], scale='tt'))
    assert_quantity_allclose(sun_B0[0], -6.54*u.deg, atol=5e-3*u.deg)
    assert_quantity_allclose(sun_B0[1], 0.88*u.deg, atol=5e-3*u.deg)


def test_B0_jpl_horizons():
    # Validate against values from JPL Horizons, which does not apply the aberration correction
    # for observer motion.
    # https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE=OBSERVER&OBJ_DATA=NO&QUANTITIES=%2714%27&COMMAND=%22Sun%22&CENTER=%27Geocentric%27&START_TIME=%222013-01-01+TT%22&STOP_TIME=%222013-12-31%22&STEP_SIZE=%221d%22
    jpl_values = {'2013-01-01': -3.034687,
                  '2013-02-01': -6.031910,
                  '2013-03-01': -7.219984,
                  '2013-04-01': -6.543092,
                  '2013-05-01': -4.168016,
                  '2013-06-01': -0.663804,
                  '2013-07-01': +2.873030,
                  '2013-08-01': +5.784693,
                  '2013-09-01': +7.194559,
                  '2013-10-01': +6.719875,
                  '2013-11-01': +4.378979,
                  '2013-12-01': +0.881068}
    sun_B0 = sun.B0(Time(list(jpl_values.keys()), scale='tt'))
    assert_quantity_allclose(sun_B0, list(jpl_values.values()) * u.deg, atol=0.01*u.arcsec, rtol=0)


def test_B0_sunspice():
    # Validate against values from SunSPICE (including calling CSPICE functions)
    # Without the aberration correction for observer motion (specify 'LT')
    #
    # IDL> load_sunspice_gen
    # IDL> cspice_str2et, '2013-01-01', et
    # IDL> cspice_subpnt, 'Intercept/Ellipsoid', 'Sun', et, 'IAU_Sun', 'LT', 'Earth', spoint, trgepc, srfvec
    # IDL> print, spclat * cspice_dpr()
    #       -3.0347784
    jpl_values = {'2013-01-01': -3.0347784,
                  '2013-02-01': -6.0319658,
                  '2013-03-01': -7.2199937,
                  '2013-04-01': -6.5430498,
                  '2013-05-01': -4.1679382,
                  '2013-06-01': -0.66371004,
                  '2013-07-01': +2.8731155,
                  '2013-08-01': +5.7847503,
                  '2013-09-01': +7.1945709,
                  '2013-10-01': +6.7198381,
                  '2013-11-01': +4.3789008,
                  '2013-12-01': +0.88096965}
    sun_B0 = sun.B0(Time(list(jpl_values.keys()), scale='utc'))
    assert_quantity_allclose(sun_B0, list(jpl_values.values()) * u.deg, atol=0.005*u.arcsec, rtol=0)


def test_L0_defaults():
    # Confirm that the output is the same when the default settings are explicitly set
    defaults = sun.L0('1992-Oct-13')
    explicit = sun.L0('1992-Oct-13',
                      light_travel_time_correction=True,
                      aberration_correction=False,
                      nearest_point=True)
    assert defaults == explicit


def test_L0_astronomical_almanac():
    # Validate against published values from the Astronomical Almanac (2013)
    aa_values = {'2013-01-01': 326.98,
                 '2013-02-01': 278.78,
                 '2013-03-01': 270.07,
                 '2013-04-01': 221.44,
                 '2013-05-01': 185.30,
                 '2013-06-01': 135.30,
                 '2013-07-01': 98.22,
                 '2013-08-01': 48.03,
                 '2013-09-01': 358.28,
                 '2013-10-01': 322.22,
                 '2013-11-01': 273.31,
                 '2013-12-01': 237.83}
    sun_L0 = sun.L0(Time(list(aa_values.keys()), scale='tt'),
                    light_travel_time_correction=True,
                    aberration_correction=True,
                    nearest_point=False)
    assert_longitude_allclose(sun_L0, list(aa_values.values()) * u.deg, atol=5e-3*u.deg)


def test_L0_jpl_horizons():
    # Validate against values from JPL Horizons, which does not apply the aberration correction
    # for observer motion.
    # https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&TABLE_TYPE=OBSERVER&OBJ_DATA=NO&QUANTITIES=%2714%27&COMMAND=%22Sun%22&CENTER=%27Geocentric%27&START_TIME=%222013-01-01+TT%22&STOP_TIME=%222013-12-31%22&STEP_SIZE=%221d%22
    jpl_values = {'2013-01-01': 326.989970,
                  '2013-02-01': 278.789894,
                  '2013-03-01': 270.072186,
                  '2013-04-01': 221.440599,
                  '2013-05-01': 185.306476,
                  '2013-06-01': 135.303097,
                  '2013-07-01': 98.221806,
                  '2013-08-01': 48.035951,
                  '2013-09-01': 358.289921,
                  '2013-10-01': 322.226009,
                  '2013-11-01': 273.315206,
                  '2013-12-01': 237.836959}
    sun_L0 = sun.L0(Time(list(jpl_values.keys()), scale='tt'),
                    light_travel_time_correction=True,
                    aberration_correction=False,
                    nearest_point=True)
    assert_longitude_allclose(sun_L0, list(jpl_values.values()) * u.deg, atol=0.01*u.arcsec)


def test_L0_sunspice():
    # Validate against values from SunSPICE (including calling CSPICE functions)

    # With the aberration correction for observer motion (specify 'LT+S')
    #
    # IDL> load_sunspice_gen
    # IDL> cspice_str2et, '2013-01-01', et
    # IDL> cspice_subpnt, 'Intercept/Ellipsoid', 'Sun', et, 'IAU_Sun', 'LT+S', 'Earth', spoint1, trgepc1, srfvec1
    # IDL> cspice_reclat, spoint1, spcrad1, spclon1, spclat1
    # IDL> print, spclon1 * cspice_dpr()
    #       -33.025998
    values1 = {'2013-01-01': -33.025998,
               '2013-02-01': -81.226108,
               '2013-03-01': -89.943817,
               '2013-04-01': -138.57536,
               '2013-05-01': -174.70941,
               '2013-06-01': +135.28726,
               '2013-07-01': +98.205970,
               '2013-08-01': +48.020071,
               '2013-09-01': -1.7260099,
               '2013-10-01': -37.789945,
               '2013-11-01': -86.700744,
               '2013-12-01': -122.17899}
    sun_L0 = sun.L0(Time(list(values1.keys()), scale='utc'),
                    light_travel_time_correction=True,
                    aberration_correction=True,
                    nearest_point=True)
    assert_longitude_allclose(sun_L0, list(values1.values()) * u.deg, atol=0.3*u.arcsec)

    # Without the aberration correction for observer motion (specify 'LT')
    #
    # IDL> cspice_subpnt, 'Intercept/Ellipsoid', 'Sun', et, 'IAU_Sun', 'LT', 'Earth', spoint2, trgepc2, srfvec2
    # IDL> cspice_reclat, spoint2, spcrad2, spclon2, spclat2
    # IDL> print, spclon2 * cspice_dpr()
    #       -33.020271
    values2 = {'2013-01-01': -33.020271,
               '2013-02-01': -81.220344,
               '2013-03-01': -89.938057,
               '2013-04-01': -138.56966,
               '2013-05-01': -174.70380,
               # '2013-06-01':  135.29281,  # skipping due to low precision in comparison
               '2013-07-01': +98.211514,
               '2013-08-01': +48.025667,
               '2013-09-01': -1.7203500,
               '2013-10-01': -37.784252,
               '2013-11-01': -86.695047,
               '2013-12-01': -122.17329}
    sun_L0 = sun.L0(Time(list(values2.keys()), scale='utc'),
                    light_travel_time_correction=True,
                    aberration_correction=False,
                    nearest_point=True)
    assert_longitude_allclose(sun_L0, list(values2.values()) * u.deg, atol=0.01*u.arcsec)

    # Without any corrections (do a straight conversion from 'HEQ' to 'Carrington')
    #
    # IDL> coord = [1.d, 0.d, 10.d]
    # IDL> convert_sunspice_lonlat, '2013-01-01', coord, 'HEQ', 'Carrington', /au, /degrees
    # IDL> print, coord
    #        1.0000000       326.89956       10.000000
    values3 = {'2013-01-01': 326.89956,
               '2013-02-01': 278.69932,
               '2013-03-01': 269.98115,
               '2013-04-01': 221.34886,
               '2013-05-01': 185.21404,
               '2013-06-01': 135.21012,
               '2013-07-01': 98.128607,
               '2013-08-01': 47.942897,
               '2013-09-01': 358.19735,
               '2013-10-01': 322.13410,
               '2013-11-01': 273.22402,
               '2013-12-01': 237.74631}
    sun_L0 = sun.L0(Time(list(values3.keys()), scale='utc'),
                    light_travel_time_correction=False,
                    aberration_correction=False,
                    nearest_point=True)
    assert_longitude_allclose(sun_L0, list(values3.values()) * u.deg, atol=0.02*u.arcsec)


def test_P():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(sun.P('1992-Oct-13'), 26.27*u.deg, atol=5e-3*u.deg)


def test_P_array_time():
    # Validate against published values from the Astronomical Almanac (2013)
    sun_P = sun.P(Time(['2013-01-01',
                        '2013-02-01',
                        '2013-03-01',
                        '2013-04-01',
                        '2013-05-01',
                        '2013-06-01',
                        '2013-07-01',
                        '2013-08-01',
                        '2013-09-01',
                        '2013-10-01',
                        '2013-11-01',
                        '2013-12-01'], scale='tt'))
    assert_quantity_allclose(sun_P, [1.98,
                                     -12.23,
                                     -21.55,
                                     -26.15,
                                     -24.11,
                                     -15.39,
                                     -2.64,
                                     10.85,
                                     21.08,
                                     25.97,
                                     24.47,
                                     16.05]*u.deg, atol=5e-3*u.deg)


def test_earth_distance():
    # Validate against a published value from Astronomical Algorithms (Meeus 1998, p.191)
    assert_quantity_allclose(sun.earth_distance('1992-Oct-13'), 0.997608*u.AU, atol=5e-7*u.AU)


def test_earth_distance_array_time():
    # Validate against published values from the Astronomical Almanac (2013)
    sunearth_distance = sun.earth_distance(Time(['2013-04-01', '2013-12-01'], scale='tt'))
    assert_quantity_allclose(sunearth_distance[0], 0.9992311*u.AU, rtol=0, atol=5e-8*u.AU)
    assert_quantity_allclose(sunearth_distance[1], 0.9861362*u.AU, rtol=0, atol=5e-8*u.AU)


def test_orientation():
    # Not currently aware of a published value to check against, so just self-check for now

    # Check the Northern Hemisphere
    angle = sun.orientation(EarthLocation(lat=40*u.deg, lon=-75*u.deg), '2017-07-18 12:00')
    assert_quantity_allclose(angle, -59.4*u.deg, atol=0.1*u.deg)

    # Check the Southern Hemisphere
    angle = sun.orientation(EarthLocation(lat=-40*u.deg, lon=-75*u.deg), '2017-02-18 13:00')
    assert_quantity_allclose(angle, -110.8*u.deg, atol=0.1*u.deg)


# Validate against published values from the Astronomical Almanac (2013, C4)
@pytest.mark.parametrize(("date", "day_fraction", "rotation_number"),
                         [('2012-12-29', 0.49, 2132),
                          ('2013-01-25', 0.83, 2133),
                          ('2013-02-22', 0.17, 2134),
                          ('2013-03-21', 0.49, 2135),
                          ('2013-04-17', 0.78, 2136),
                          ('2013-05-15', 0.02, 2137),
                          ('2013-06-11', 0.22, 2138),
                          ('2013-07-08', 0.42, 2139),
                          ('2013-08-04', 0.63, 2140),
                          ('2013-08-31', 0.87, 2141),
                          ('2013-09-28', 0.14, 2142),
                          ('2013-10-25', 0.43, 2143),
                          ('2013-11-21', 0.73, 2144),
                          ('2013-12-19', 0.05, 2145),
                          ('2014-01-15', 0.38, 2146),
                          ('2014-02-11', 0.73, 2147)])
def test_carrington_rotation_number(date, day_fraction, rotation_number):
    assert_allclose(sun.carrington_rotation_number(Time(date, scale='tt') + day_fraction*u.day),
                    rotation_number, rtol=0, atol=2e-4)


@pytest.mark.parametrize(("crot", "julian_days"),
                         [(1, 2398167.4),
                          (2, 2398194.6756),
                          (1860, 2448872.027509),
                          ])
def test_carrington_rotation_starttime(crot, julian_days):
    # Stated precision in the docstring is 0.11 seconds
    with warnings.catch_warnings():
        # Filter warnings caused by very old dates
        warnings.filterwarnings("ignore", category=ErfaWarning)
        assert_quantity_allclose(sun.carrington_rotation_time(crot).tt.jd * u.day,
                                 julian_days * u.day, atol=0.11*u.s)


@pytest.mark.parametrize(("crot", "longitude", "crot_fractional"),
                         [(2000, 360, 2000.0),
                          (2000.0, 270, 2000.25)
                          ])
def test_carrington_rotation_time_longitude(crot, longitude, crot_fractional):
    assert sun.carrington_rotation_time(crot*u.one, longitude*u.deg) == \
        sun.carrington_rotation_time(crot_fractional*u.one)


@pytest.mark.parametrize(("crot", "longitude", "crot_fractional"),
                         [
                             (np.array([2000, 2000]), np.array(
                                 [180, 90]), np.array([2000.5, 2000.75])),
                             (2000, np.array([180, 90]), np.array([2000.5, 2000.75])),
                             (np.array([2000, 2000]), 180, np.array([2000.5, 2000.5]))
])
def test_carrington_rotation_time_longitude_numpy(crot, longitude, crot_fractional):
    assert all(sun.carrington_rotation_time(crot*u.one, longitude*u.deg) ==
               sun.carrington_rotation_time(crot_fractional*u.one))


@pytest.mark.parametrize(("crot", "longitude", "expected_error"),
    [
        (2000, 0, "Carrington longitude(s) must be > 0 degrees and <= 360 degrees."),
        (2000, -10, "Carrington longitude(s) must be > 0 degrees and <= 360 degrees."),
        (2000, 400, "Carrington longitude(s) must be > 0 degrees and <= 360 degrees."),
        (2000.5, 180, "Carrington rotation number(s) must be integral if `longitude` is provided.")
    ]
)
def test_carrington_rotation_time_longitude_err(crot, longitude, expected_error):
    with pytest.raises(ValueError, match=re.escape(expected_error)):
        sun.carrington_rotation_time(crot*u.one, longitude*u.deg)


def test_carrington_rotation_roundtrip():
    t = Time('2010-1-1')
    crot = sun.carrington_rotation_number(t)
    t_roundtrip = sun.carrington_rotation_time(crot)
    dt = t - t_roundtrip
    # Stated precision in the docstring is 0.11 seconds
    assert_quantity_allclose(dt.to(u.s), 0*u.s, atol=0.11*u.s)


def test_carrington_rotation_str():
    # Check that by default a human parseable string is returned
    t = sun.carrington_rotation_time(2210)
    assert str(t) == '2018-10-26 20:48:16.137'


# For 2024 Apr 8, at 29.6 deg N, 98.5 deg W, one eclipse calculator has:
#   Partial eclipse begins: 17:14:49 UTC
#   Totality begins: 18:33:44 UTC
#   Maximum eclipse: 18:34:38 UTC
#   Totality ends: 18:35:32 UTC
#   Partial eclipse ends: 19:56:00 UTC
# https://www.timeanddate.com/eclipse/in/@29.6,-98.5?iso=20240408

@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore:Tried to get polar motions for times after IERS data is valid.")
@pytest.mark.filterwarnings("ignore:.*times are outside of range covered by IERS table.")
def test_eclipse_amount(use_DE440s):
    location = EarthLocation.from_geodetic(-98.5*u.deg, 29.6*u.deg)

    # We use the mean lunar radius (the default) for penumbral contacts
    assert sun.eclipse_amount(location.get_itrs(Time('2024-04-08 17:14:48'))) == 0
    assert sun.eclipse_amount(location.get_itrs(Time('2024-04-08 17:14:53'))) > 0
    assert sun.eclipse_amount(location.get_itrs(Time('2024-04-08 19:55:58'))) > 0
    assert sun.eclipse_amount(location.get_itrs(Time('2024-04-08 19:56:01'))) == 0


@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore:Tried to get polar motions for times after IERS data is valid.")
@pytest.mark.filterwarnings("ignore:.*times are outside of range covered by IERS table.")
def test_eclipse_amount_minimum(use_DE440s):
    location = EarthLocation.from_geodetic(-98.5*u.deg, 29.6*u.deg)

    # We use the mean minimum lunar radius for umbral contacts
    assert sun.eclipse_amount(location.get_itrs(Time('2024-04-08 18:33:43')), moon_radius="minimum") < 1
    assert sun.eclipse_amount(location.get_itrs(Time('2024-04-08 18:33:45')), moon_radius="minimum") == 1
    assert sun.eclipse_amount(location.get_itrs(Time('2024-04-08 18:35:31')), moon_radius="minimum") == 1
    assert sun.eclipse_amount(location.get_itrs(Time('2024-04-08 18:35:35')), moon_radius="minimum") < 1
