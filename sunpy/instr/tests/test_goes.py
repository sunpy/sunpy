from __future__ import absolute_import

import copy
import datetime
import pytest

import numpy as np
from astropy.units.quantity import Quantity
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_array_equal, assert_almost_equal
from pandas.util.testing import assert_frame_equal
import astropy.units as u

from sunpy.time import TimeRange
from sunpy import lightcurve
from sunpy.instr import goes

# Define input variables to be used in test functions for
# _goes_chianti_tem.
LONGFLUX = Quantity([7e-6], unit="W/m**2")
SHORTFLUX = Quantity([7e-7], unit="W/m**2")
DATE = "2014-04-16"

@pytest.mark.remote_data
def test_goes_event_list():
    # Set a time range to search
    trange = TimeRange('2011-06-07 00:00', '2011-06-08 00:00')
    # Test case where GOES class filter is applied
    result = goes.get_goes_event_list(trange, goes_class_filter='M1')
    assert type(result) == list
    assert type(result[0]) == dict
    assert type(result[0]['event_date'] == str)
    assert type(result[0]['goes_location'] == tuple)
    assert type(result[0]['peak_time'] == datetime.datetime)
    assert type(result[0]['start_time'] == datetime.datetime)
    assert type(result[0]['end_time'] == datetime.datetime)
    assert type(result[0]['goes_class'] == str)
    assert type(result[0]['noaa_active_region'] == int)
    assert result[0]['event_date'] == '2011-06-07'
    assert result[0]['goes_location'] == (54, -21)
    assert result[0]['start_time'] == datetime.datetime(2011, 6, 7, 6, 16)
    assert result[0]['peak_time'] == datetime.datetime(2011, 6, 7, 6, 41)
    assert result[0]['end_time'] == datetime.datetime(2011, 6, 7, 6, 59)
    assert result[0]['goes_class'] == 'M2.5'
    assert result[0]['noaa_active_region'] == 11226
    # Test case where GOES class filter not applied
    result = goes.get_goes_event_list(trange)
    assert type(result) == list
    assert type(result[0]) == dict
    assert type(result[0]['event_date'] == str)
    assert type(result[0]['goes_location'] == tuple)
    assert type(result[0]['peak_time'] == datetime.datetime)
    assert type(result[0]['start_time'] == datetime.datetime)
    assert type(result[0]['end_time'] == datetime.datetime)
    assert type(result[0]['goes_class'] == str)
    assert type(result[0]['noaa_active_region'] == int)
    assert result[0]['event_date'] == '2011-06-07'
    assert result[0]['goes_location'] == (54, -21)
    assert result[0]['start_time'] == datetime.datetime(2011, 6, 7, 6, 16)
    assert result[0]['peak_time'] == datetime.datetime(2011, 6, 7, 6, 41)
    assert result[0]['end_time'] == datetime.datetime(2011, 6, 7, 6, 59)
    assert result[0]['goes_class'] == 'M2.5'
    assert result[0]['noaa_active_region'] == 11226


@pytest.mark.remote_data
def test_calculate_temperature_em():
    # Create GOESLightcurve object, then create new one with
    # temperature & EM using with calculate_temperature_em().
    goeslc = lightcurve.GOESLightCurve.create("2014-01-01 00:00", "2014-01-01 01:00")
    goeslc_new = goes.calculate_temperature_em(goeslc)
    # Test correct exception is raised if a GOESLightCurve object is
    # not inputted.
    with pytest.raises(TypeError):
        goes.calculate_temperature_em([])
    # Find temperature and EM manually with _goes_chianti_tem()
    temp, em = goes._goes_chianti_tem(
        Quantity(goeslc.data.xrsb, unit='W/m**2'),
        Quantity(goeslc.data.xrsa, unit='W/m**2'),
        satellite=int(goeslc.meta["TELESCOP"].split()[1]), date="2014-01-01")
    # Check that temperature and EM arrays from _goes_chianti_tem()
    # are same as those in new GOESLightcurve object.
    assert goeslc_new.data.temperature.all() == temp.value.all()
    assert goeslc_new.data.em.all() == em.value.all()
    # Check rest of data frame of new GOESLightCurve object is same
    # as that in original object.
    goeslc_revert = copy.deepcopy(goeslc_new)
    del goeslc_revert.data["temperature"]
    del goeslc_revert.data["em"]
    assert_frame_equal(goeslc_revert.data, goeslc.data)

@pytest.mark.remote_data
def test_goes_chianti_tem_errors():
    # Define input variables.
    ratio = SHORTFLUX/LONGFLUX
    shortflux_toomany = Quantity(
        np.append(SHORTFLUX.value, SHORTFLUX.value[0]), unit="W/m**2")
    shortflux_toosmall = copy.deepcopy(SHORTFLUX)
    shortflux_toosmall.value[0] = -1
    shortflux_toobig = copy.deepcopy(SHORTFLUX)
    shortflux_toobig.value[0] = 1
    temp_test = Quantity(np.zeros(len(LONGFLUX))+10, unit="MK")
    temp_test_toomany = Quantity(np.append(temp_test.value, 0), unit="MK")
    temp_test_toosmall = copy.deepcopy(temp_test)
    temp_test_toosmall.value[0] = -1
    temp_test_toobig = copy.deepcopy(temp_test)
    temp_test_toobig.value[0] = 101
    # First test correct exceptions are raised if incorrect inputs are
    # entered.
    with pytest.raises(ValueError):
        temp, em = goes._goes_chianti_tem(LONGFLUX, SHORTFLUX, satellite=-1)
    with pytest.raises(ValueError):
        temp, em = goes._goes_chianti_tem(LONGFLUX, shortflux_toomany)
    with pytest.raises(ValueError):
        temp = goes._goes_get_chianti_temp(ratio, satellite=-1)
    with pytest.raises(ValueError):
        temp, em = goes._goes_chianti_tem(LONGFLUX, SHORTFLUX,
                                         abundances="Neither")
    with pytest.raises(ValueError):
        temp = goes._goes_get_chianti_temp(ratio, abundances="Neither")
    with pytest.raises(ValueError):
        temp, em = goes._goes_chianti_tem(LONGFLUX, shortflux_toobig)
    with pytest.raises(ValueError):
        em = goes._goes_get_chianti_em(LONGFLUX, temp_test, satellite=-1)
    with pytest.raises(ValueError):
        em = goes._goes_get_chianti_em(LONGFLUX, temp_test,
                                       abundances="Neither")
    with pytest.raises(ValueError):
        em = goes._goes_get_chianti_em(LONGFLUX, temp_test,
                                       abundances="Neither")
    with pytest.raises(ValueError):
        em = goes._goes_get_chianti_em(LONGFLUX, temp_test_toosmall)
    with pytest.raises(ValueError):
        em = goes._goes_get_chianti_em(LONGFLUX, temp_test_toobig)

@pytest.mark.remote_data
def test_goes_chianti_tem_case1():
    # test case 1: satellite > 7, abundances = coronal
    temp1, em1 = goes._goes_chianti_tem(LONGFLUX, SHORTFLUX, satellite=15,
                                       date=DATE)
    np.testing.assert_allclose(temp1, Quantity([11.28], unit="MK"), rtol=0.01)
    assert all(em1 < Quantity([4.79e+48], unit="1/cm**3")) and \
      em1 > Quantity([4.78e+48], unit="1/cm**3")

@pytest.mark.remote_data
def test_goes_chianti_tem_case2():
    # test case 2: satellite > 7, abundances = photospheric
    temp2, em2 = goes._goes_chianti_tem(LONGFLUX, SHORTFLUX, satellite=15,
                                       date=DATE, abundances="photospheric")
    assert all(temp2 < Quantity([10.25], unit="MK")) and \
      all(temp2 > Quantity([10.24], unit="MK"))
    assert all(em2 < Quantity([1.12e+49], unit="1/cm**3")) and \
      all(em2 > Quantity([1.11e+49], unit="1/cm**3"))

@pytest.mark.remote_data
def test_goes_chianti_tem_case3():
    # test case 3: satellite < 8 and != 6, abundances = coronal
    temp3, em3 = goes._goes_chianti_tem(LONGFLUX, SHORTFLUX, satellite=5,
                                       date=DATE,
                                       abundances="coronal")
    assert all(temp3 < Quantity([11.43], unit="MK")) and \
      all(temp3 > Quantity([11.42], unit="MK"))
    assert all(em3 < Quantity([3.85e+48], unit="1/cm**3")) and \
      all(em3 > Quantity([3.84e+48], unit="1/cm**3"))

@pytest.mark.remote_data
def test_goes_chianti_tem_case4():
    # test case 4: satellite < 8 and != 6, abundances = photospheric
    temp4, em4 = goes._goes_chianti_tem(LONGFLUX, SHORTFLUX, satellite=5,
                                       date=DATE,
                                       abundances="photospheric")
    assert all(temp4 < Quantity([10.42], unit="MK")) and \
      all(temp4 > Quantity([10.41], unit="MK"))
    assert all(em4 < Quantity(8.81e+48, unit="1/cm**3")) and \
      all(em4 > Quantity(8.80e+48, unit="1/cm**3"))

@pytest.mark.remote_data
def test_goes_chianti_tem_case5():
    # test case 5: satellite = 6, date < 1983-06-28, abundances = coronal
    temp5, em5 = goes._goes_chianti_tem(LONGFLUX, SHORTFLUX, satellite=6,
                                       date="1983-06-27",
                                       abundances="coronal")
    assert all(temp5 < Quantity(12.30, unit="MK")) and \
      all(temp5 > Quantity(12.29, unit="MK"))
    assert all(em5 < Quantity(3.13e+48, unit="1/cm**3")) and \
      all(em5 > Quantity(3.12e+48, unit="1/cm**3"))

@pytest.mark.remote_data
def test_goes_chianti_tem_case6():
    # test case 6: satellite = 6, date < 1983-06-28, abundances = photospheric
    temp6, em6 = goes._goes_chianti_tem(LONGFLUX, SHORTFLUX, satellite=6,
                                       date="1983-06-27",
                                       abundances="photospheric")
    assert all(temp6 < Quantity(11.44, unit="MK")) and \
      all(temp6 > Quantity(11.43, unit="MK"))
    assert all(em6 < Quantity(6.74e+48, unit="1/cm**3")) and \
      all(em6 > Quantity(6.73e+48, unit="1/cm**3"))

@pytest.mark.remote_data
def test_goes_chianti_tem_case7():
    # test case 7: satellite = 6, date > 1983-06-28, abundances = coronal
    temp7, em7 = goes._goes_chianti_tem(LONGFLUX, SHORTFLUX, satellite=6,
                                       date=DATE,
                                       abundances="coronal")
    assert all(temp7 < Quantity(11.34, unit="MK")) and \
      all(temp7 > Quantity(11.33, unit="MK"))
    assert all(em7 < Quantity(4.08e+48, unit="1/cm**3")) and \
      all(em7 > Quantity(4.07e+48, unit="1/cm**3"))

@pytest.mark.remote_data
def test_goes_chianti_tem_case8():
    # test case 8: satellite = 6, date > 1983-06-28, abundances = photospheric
    temp8, em8 = goes._goes_chianti_tem(LONGFLUX, SHORTFLUX, satellite=6,
                                       date=DATE,
                                       abundances="photospheric")
    assert all(temp8 < Quantity(10.36, unit="MK")) and \
      all(temp8 > Quantity(10.35, unit="MK"))
    assert all(em8 < Quantity(9.39e+48, unit="1/cm**3")) and \
      all(em8 > Quantity(9.38e+48, unit="1/cm**3"))

@pytest.mark.remote_data
def test_calculate_radiative_loss_rate():
    # Define input variables.
    goeslc_input = lightcurve.GOESLightCurve.create("2014-01-01 00:00:00",
                                            "2014-01-01 00:00:10")
    not_goeslc = []
    goeslc_no_em = goes.calculate_temperature_em(goeslc_input)
    del goeslc_no_em.data["em"]

    # Check correct exceptions are raised to incorrect inputs
    with pytest.raises(TypeError):
        goes_test = goes.calculate_radiative_loss_rate(not_goeslc)

    # Check function gives correct results.
    # Test case 1: GOESLightCurve object with only flux data
    goeslc_test = goes.calculate_radiative_loss_rate(goeslc_input)
    goeslc_expected = goes.calculate_temperature_em(goeslc_input)
    goeslc_expected.data["rad_loss_rate"] = \
      np.array([5.44914366e+19, 5.44914366e+19, 5.43465905e+19,
                5.38282295e+19, 5.42019309e+19])
    assert_frame_equal(goeslc_test.data, goeslc_expected.data)

    # Test case 2: GOESLightCurve object with flux and temperature
    # data, but no EM data.
    goes_test = goes.calculate_radiative_loss_rate(goeslc_no_em)
    assert_frame_equal(goeslc_test.data, goeslc_expected.data)

@pytest.mark.remote_data
def test_calc_rad_loss_errors():
    # Define input variables
    temp = 11.0 * Quantity(np.ones(6), unit="MK")
    em = 4.0e+48 * Quantity(np.ones(6), unit="1/cm**3")
    obstime = np.array([datetime.datetime(2014, 1, 1, 0, 0, 0),
                        datetime.datetime(2014, 1, 1, 0, 0, 2),
                        datetime.datetime(2014, 1, 1, 0, 0, 4),
                        datetime.datetime(2014, 1, 1, 0, 0, 6),
                        datetime.datetime(2014, 1, 1, 0, 0, 8),
                        datetime.datetime(2014, 1, 1, 0, 0, 10)], dtype=object)
    temp_toolong = Quantity(np.append(temp.value, 0), unit="MK")
    obstime_toolong =  np.array([datetime.datetime(2014, 1, 1, 0, 0, 0),
                        datetime.datetime(2014, 1, 1, 0, 0, 2),
                        datetime.datetime(2014, 1, 1, 0, 0, 4),
                        datetime.datetime(2014, 1, 1, 0, 0, 6),
                        datetime.datetime(2014, 1, 1, 0, 0, 8),
                        datetime.datetime(2014, 1, 1, 0, 0, 10),
                        datetime.datetime(2014, 1, 1, 0, 0, 12)], dtype=object)
    obstime_nonchrono = copy.deepcopy(obstime)
    obstime_nonchrono[1] = obstime[-1]
    obstime_nonchrono[-1] = obstime[1]
    obstime_notdatetime = copy.deepcopy(obstime)
    obstime_notdatetime[0] = 1
    temp_outofrange = Quantity([101, 11.0, 11.0, 11.0, 11.0, 11.0], unit="MK")
    # Ensure correct exceptions are raised.
    with pytest.raises(ValueError):
        rad_loss_test = goes._calc_rad_loss(temp_toolong, em, obstime)
    with pytest.raises(ValueError):
        rad_loss_test = goes._calc_rad_loss(temp_outofrange, em, obstime)
    with pytest.raises(IOError):
        rad_loss_test = goes._calc_rad_loss(temp, em, obstime_toolong)
    with pytest.raises(TypeError):
        lx_test = goes._calc_rad_loss(temp, em, obstime_notdatetime)
    with pytest.raises(ValueError):
        rad_loss_test = goes._calc_rad_loss(temp, em, obstime_nonchrono)

@pytest.mark.remote_data
def test_calc_rad_loss_nokwags():
    # Define input variables
    temp = Quantity([11.0, 11.0, 11.0, 11.0, 11.0, 11.0], unit="MK")
    em = Quantity([4.0e+48, 4.0e+48, 4.0e+48, 4.0e+48, 4.0e+48, 4.0e+48],
                  unit="1/cm**3")
    obstime = np.array([datetime.datetime(2014, 1, 1, 0, 0, 0),
                        datetime.datetime(2014, 1, 1, 0, 0, 2),
                        datetime.datetime(2014, 1, 1, 0, 0, 4),
                        datetime.datetime(2014, 1, 1, 0, 0, 6),
                        datetime.datetime(2014, 1, 1, 0, 0, 8),
                        datetime.datetime(2014, 1, 1, 0, 0, 10)], dtype=object)
    # Test output is correct when no kwags are set.
    rad_loss_test = goes._calc_rad_loss(temp[:2], em[:2])
    rad_loss_expected = {"rad_loss_rate":
                         3.01851392e+19 * Quantity(np.ones(2), unit="J/s")}
    assert sorted(rad_loss_test.keys()) == sorted(rad_loss_expected.keys())
    assert_quantity_allclose(rad_loss_test["rad_loss_rate"],
                       rad_loss_expected["rad_loss_rate"], rtol=0.01)

@pytest.mark.remote_data
def test_calc_rad_loss_obstime():
    # Define input variables
    temp = Quantity([11.0, 11.0, 11.0, 11.0, 11.0, 11.0], unit="MK")
    em = Quantity([4.0e+48, 4.0e+48, 4.0e+48, 4.0e+48, 4.0e+48, 4.0e+48],
                  unit="1/cm**3")
    obstime = np.array([datetime.datetime(2014, 1, 1, 0, 0, 0),
                        datetime.datetime(2014, 1, 1, 0, 0, 2),
                        datetime.datetime(2014, 1, 1, 0, 0, 4),
                        datetime.datetime(2014, 1, 1, 0, 0, 6),
                        datetime.datetime(2014, 1, 1, 0, 0, 8),
                        datetime.datetime(2014, 1, 1, 0, 0, 10)], dtype=object)
    # Test output is correct when obstime and cumulative kwargs are set.
    rad_loss_test = goes._calc_rad_loss(temp, em, obstime)
    rad_loss_expected = {
        "rad_loss_rate": 3.01851392e+19 * Quantity(np.ones(6), unit="J/s"),
        "rad_loss_int": Quantity(3.01851392e+20, unit="J"),
        "rad_loss_cumul": Quantity([6.03702783e+19, 1.20740557e+20,
                                    1.81110835e+20, 2.41481113e+20,
                                    3.01851392e+20], unit="J")
        }
    assert sorted(rad_loss_test.keys()) == sorted(rad_loss_expected.keys())
    assert_quantity_allclose(rad_loss_test["rad_loss_rate"],
                       rad_loss_expected["rad_loss_rate"], rtol=0.0001)
    assert_quantity_allclose(rad_loss_test["rad_loss_int"],
                       rad_loss_expected["rad_loss_int"], rtol=0.0001)
    assert_quantity_allclose(rad_loss_test["rad_loss_cumul"],
                       rad_loss_expected["rad_loss_cumul"], rtol=0.0001)

@pytest.mark.remote_data
def test_calculate_xray_luminosity():
    # Check correct exceptions are raised to incorrect inputs
    not_goeslc = []
    with pytest.raises(TypeError):
        goes_test = goes.calculate_xray_luminosity(not_goeslc)
    # Check function gives correct results.
    goeslc_input = lightcurve.GOESLightCurve.create("2014-01-01 00:00:00",
                                            "2014-01-01 00:00:10")
    goeslc_test = goes.calculate_xray_luminosity(goeslc_input)
    goeslc_expected = copy.deepcopy(goeslc_input)
    goeslc_expected.data["luminosity_xrsa"] = \
      Quantity(np.array([2.49831950e+16, 2.49831950e+16, 2.49831950e+16,
                         2.52864004e+16, 2.49831950e+16], dtype="float32"),
                         unit="J/s")
    goeslc_expected.data["luminosity_xrsb"] = \
      Quantity(np.array([9.54399250e+17, 9.54399250e+17, 9.52985195e+17,
                         9.52985195e+17, 9.51571139e+17], dtype="float32"),
                         unit="J/s")
    assert_frame_equal(goeslc_test.data, goeslc_expected.data,
                       check_less_precise=True)

def test_goes_lx_errors():
    # Define input values of flux and time.
    longflux = 7e-6 * Quantity(np.ones(6), unit="W/m**2")
    shortflux = 7e-7 * Quantity(np.ones(6), unit="W/m**2")
    obstime = np.array([datetime.datetime(2014, 1, 1, 0, 0, 0),
                        datetime.datetime(2014, 1, 1, 0, 0, 2),
                        datetime.datetime(2014, 1, 1, 0, 0, 4),
                        datetime.datetime(2014, 1, 1, 0, 0, 6),
                        datetime.datetime(2014, 1, 1, 0, 0, 8),
                        datetime.datetime(2014, 1, 1, 0, 0, 10)], dtype=object)
    longflux_toolong = Quantity(np.append(longflux.value, 0), unit=longflux.unit)
    obstime_nonchrono = copy.deepcopy(obstime)
    obstime_nonchrono[1] = obstime[-1]
    obstime_notdatetime = copy.deepcopy(obstime)
    obstime_notdatetime[0] = 1
    # Ensure correct exceptions are raised.
    with pytest.raises(ValueError):
        lx_test = goes._goes_lx(longflux_toolong, shortflux, obstime)
    with pytest.raises(TypeError):
        lx_test = goes._goes_lx(longflux, shortflux, obstime_notdatetime)
    with pytest.raises(ValueError):
        lx_test = goes._goes_lx(longflux, shortflux, obstime_nonchrono)

def test_goes_lx_nokwargs():
    # Define input values of flux and time.
    longflux = Quantity([7e-6, 7e-6, 7e-6, 7e-6, 7e-6, 7e-6], unit="W/m**2")
    shortflux = Quantity([7e-7, 7e-7, 7e-7, 7e-7, 7e-7, 7e-7], unit="W/m**2")
    # Test output when no kwargs are set.
    lx_test = goes._goes_lx(longflux[:2], shortflux[:2])
    lx_expected = {"longlum": Quantity([1.98649103e+18, 1.98649103e+18],
                                       unit="W"),
                   "shortlum": Quantity([1.98649103e+17, 1.98649103e+17],
                                        unit="W")}
    assert sorted(lx_test.keys()) == sorted(lx_expected.keys())
    assert_quantity_allclose(lx_test["longlum"], lx_expected["longlum"], rtol=0.1)
    assert_quantity_allclose(lx_test["shortlum"], lx_expected["shortlum"],
                       rtol=0.1)

def test_goes_lx_date():
    # Define input values of flux and time.
    longflux = Quantity([7e-6, 7e-6, 7e-6, 7e-6, 7e-6, 7e-6], unit="W/m**2")
    shortflux = Quantity([7e-7, 7e-7, 7e-7, 7e-7, 7e-7, 7e-7], unit="W/m**2")
    # Test output when date kwarg is set.
    lx_test = goes._goes_lx(longflux[:2], shortflux[:2], date="2014-04-21")
    lx_expected = {"longlum": Quantity([1.98649103e+18, 1.98649103e+18],
                                       unit="W"),
                   "shortlum": Quantity([1.98649103e+17, 1.98649103e+17],
                                        unit="W")}
    assert sorted(lx_test.keys()) == sorted(lx_expected.keys())
    assert_quantity_allclose(lx_test["longlum"], lx_expected["longlum"], rtol=0.001)
    assert_quantity_allclose(lx_test["shortlum"], lx_expected["shortlum"],
                       rtol=0.001)

def test_goes_lx_obstime():
    # Define input values of flux and time.
    longflux = Quantity([7e-6, 7e-6, 7e-6, 7e-6, 7e-6, 7e-6], unit="W/m**2")
    shortflux = Quantity([7e-7, 7e-7, 7e-7, 7e-7, 7e-7, 7e-7], unit="W/m**2")
    obstime = np.array([datetime.datetime(2014, 1, 1, 0, 0, 0),
                        datetime.datetime(2014, 1, 1, 0, 0, 2),
                        datetime.datetime(2014, 1, 1, 0, 0, 4),
                        datetime.datetime(2014, 1, 1, 0, 0, 6),
                        datetime.datetime(2014, 1, 1, 0, 0, 8),
                        datetime.datetime(2014, 1, 1, 0, 0, 10)], dtype=object)
    # Test output when obstime and cumulative kwargs are set.
    lx_test = goes._goes_lx(longflux, shortflux, obstime)
    lx_expected = {
        "longlum": 1.96860565e+18 * Quantity(np.ones(6), unit='W'),
        "shortlum": 1.96860565e+17 * Quantity(np.ones(6), unit='W'),
        "longlum_int": Quantity([1.96860565e+19], unit="J"),
        "shortlum_int": Quantity([1.96860565e+18], unit="J"),
        "longlum_cumul": Quantity([3.93721131e+18, 7.87442262e+18,
                                   1.18116339e+19, 1.57488452e+19,
                                   1.96860565e+19], unit="J"),
        "shortlum_cumul": Quantity([3.93721131e+17, 7.87442262e+17,
                                    1.18116339e+18, 1.57488452e+18,
                                    1.96860565e+18], unit="J")}
    assert sorted(lx_test.keys()) == sorted(lx_expected.keys())
    assert_quantity_allclose(lx_test["longlum"], lx_expected["longlum"], rtol=0.1)
    assert_quantity_allclose(lx_test["shortlum"], lx_expected["shortlum"],
                       rtol=0.1)
    assert_quantity_allclose(lx_test["longlum_int"], lx_expected["longlum_int"],
                       rtol=0.1)
    assert_quantity_allclose(lx_test["shortlum_int"], lx_expected["shortlum_int"],
                       rtol=0.1)
    assert_quantity_allclose(lx_test["longlum_cumul"], lx_expected["longlum_cumul"],
                       rtol=0.1)
    assert_quantity_allclose(lx_test["shortlum_cumul"],
                       lx_expected["shortlum_cumul"], rtol=0.1)

def test_flux_to_classletter():
    """Test converting fluxes into a class letter"""
    fluxes = Quantity(10**(-np.arange(9, 2., -1)), 'W/m**2')
    classesletter = ['A', 'A', 'B', 'C', 'M', 'X', 'X']
    calculated_classesletter = [goes.flux_to_flareclass(f)[0] for f in fluxes]
    calculated_classnumber = [float(goes.flux_to_flareclass(f)[1:]) for f in fluxes]
    assert_array_equal(classesletter, calculated_classesletter)
    assert_array_equal([0.1, 1, 1, 1, 1, 1, 10], calculated_classnumber)
    # now test the Examples
    assert goes.flux_to_flareclass(1e-08 * u.watt/u.m**2) == 'A1'
    assert goes.flux_to_flareclass(0.00682 * u.watt/u.m**2) == 'X68.2'
    assert goes.flux_to_flareclass(7.8e-09 * u.watt/u.m**2) == 'A0.78'
    assert goes.flux_to_flareclass(0.00024 * u.watt/u.m**2) == 'X2.4'
    assert goes.flux_to_flareclass(4.7e-06 * u.watt/u.m**2) == 'C4.7'
    assert goes.flux_to_flareclass(6.9e-07 * u.watt/u.m**2) == 'B6.9'
    assert goes.flux_to_flareclass(2.1e-05 * u.watt/u.m**2) == 'M2.1'

def test_class_to_flux():
    classes = ['A3.49', 'A0.23', 'M1', 'X2.3', 'M5.8', 'C2.3', 'B3.45', 'X20']
    results = Quantity([3.49e-8, 2.3e-9, 1e-5, 2.3e-4, 5.8e-5, 2.3e-6, 3.45e-7, 2e-3], 'W/m2')
    for c, r in zip(classes, results):
        assert_almost_equal(r.value, goes.flareclass_to_flux(c).value)

def test_joint_class_to_flux():
    classes = ['A3.49', 'A0.23', 'M1', 'X2.3', 'M5.8', 'C2.3', 'B3.45', 'X20']
    for c in classes:
        assert c == goes.flux_to_flareclass(goes.flareclass_to_flux(c))

# TODO add a test to check for raising error
