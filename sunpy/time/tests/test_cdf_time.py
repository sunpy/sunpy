"""
Tests for parsing the CDF epoch formats with `sunpy.time.parse_time`.

The reference values here are cross-checked in two independent ways, so that a
mistake in how the epochs are defined cannot pass unnoticed:

* Against the definitions of the epoch types given in the `CDF User's Guide
  <https://spdf.gsfc.nasa.gov/pub/software/cdf/doc/cdf_User_Guide.pdf>`__, which
  are simple enough to check by hand and are quoted in the tests below.
* Against ``cdflib.epochs.CDFepoch``, which is an implementation written from
  those definitions directly rather than on top of `astropy.time.Time`, so it is
  a genuinely separate answer and not a restatement of the code under test.
"""
import sys
import warnings
from unittest import mock

import numpy as np
import pytest

import astropy.table
import astropy.units as u
from astropy.time import Time

from sunpy.time import is_time, parse_time

# The whole module needs cdflib, which is an optional dependency.
cdflib = pytest.importorskip("cdflib")
CDFepoch = cdflib.epochs.CDFepoch

# 2024-01-01T00:00:00 UTC expressed in each of the three CDF epoch types.
# CDF_EPOCH is milliseconds since 0000-01-01, CDF_EPOCH16 is (seconds since
# 0000-01-01) + 1j * (picoseconds within that second), and CDF_TIME_TT2000 is
# nanoseconds since 2000-01-01T12:00:00 TT.
EPOCH_2024 = 63871286400000.0
EPOCH16_2024 = 63871286400 + 0j
TT2000_2024 = 757339269184000000
ISO_2024 = '2024-01-01T00:00:00.000000000'

ALL_FORMATS = [
    ('cdf_epoch', EPOCH_2024, 'utc'),
    ('cdf_epoch16', EPOCH16_2024, 'utc'),
    ('cdf_tt2000', TT2000_2024, 'tt'),
]


@pytest.mark.parametrize('value', [EPOCH_2024, EPOCH16_2024, TT2000_2024])
def test_reference_values_match_cdflib(value):
    # cdflib.epochs.CDFepoch is an independent implementation, so agreeing with it
    # confirms the constants above really are 2024-01-01T00:00:00 UTC.
    assert CDFepoch.encode(value).startswith('2024-01-01T00:00:00')


@pytest.mark.parametrize(('format', 'value', 'scale'), ALL_FORMATS)
def test_parse_time_cdf_scalar(format, value, scale):
    t = parse_time(value, format=format)
    assert isinstance(t, Time)
    assert t.isscalar
    assert t.format == format
    # The scale is fixed by the epoch type: only TT2000 is measured in TT.
    assert t.scale == scale
    # Converting to other formats must give the same instant, which is what
    # catches a definition that is internally self-consistent but simply wrong.
    assert t.utc.isot == ISO_2024
    assert t.utc.jd == 2460310.5
    assert t.utc.mjd == 60310.0
    assert t.utc.unix == 1704067200.0
    assert t.utc.datetime64 == np.datetime64('2024-01-01T00:00:00')


@pytest.mark.parametrize(('format', 'value', 'scale'), ALL_FORMATS)
def test_parse_time_cdf_roundtrips_through_format(format, value, scale):
    # Asking astropy for the value back in the same format returns what went in.
    t = parse_time(value, format=format)
    assert getattr(t, format) == pytest.approx(np.real(value))


@pytest.mark.parametrize(('format', 'value', 'scale'), ALL_FORMATS)
@pytest.mark.parametrize('container', [list, np.array])
def test_parse_time_cdf_array(format, value, scale, container):
    # Lists and arrays reach the conversion through different branches of the
    # convert_time single dispatch, so check both.
    t = parse_time(container([value, value]), format=format)
    assert t.shape == (2,)
    assert t.format == format
    assert list(t.utc.isot) == [ISO_2024, ISO_2024]


@pytest.mark.parametrize(('format', 'value', 'scale'), ALL_FORMATS)
def test_parse_time_cdf_keeps_shape(format, value, scale):
    # A CDF time variable can have more than one dimension.
    t = parse_time(np.full((2, 3), value), format=format)
    assert t.shape == (2, 3)
    assert np.all(t.utc.isot == ISO_2024)


@pytest.mark.parametrize(('format', 'value', 'scale'), ALL_FORMATS)
def test_parse_time_cdf_pandas_and_column(format, value, scale):
    pandas = pytest.importorskip("pandas")
    expected = parse_time([value, value], format=format)
    for obj in (pandas.Series([value, value]), astropy.table.Column([value, value])):
        assert np.all(parse_time(obj, format=format).utc.isot == expected.utc.isot)


@pytest.mark.parametrize(('format', 'value'), [('cdf_epoch', 0.0),
                                               ('cdf_epoch16', 0 + 0j),
                                               ('cdf_tt2000', 0)])
def test_parse_time_cdf_epoch_zero_points(format, value):
    # Both CDF_EPOCH and CDF_EPOCH16 are counted from midnight on 0000-01-01 UTC,
    # while CDF_TIME_TT2000 is counted from 2000-01-01T12:00:00 TT (i.e. J2000).
    expected = ('2000-01-01 12:00:00', 'tt') if format == 'cdf_tt2000' else ('0000-01-01 00:00:00', 'utc')
    assert parse_time(value, format=format).jd == Time(expected[0], scale=expected[1], format='iso').jd


@pytest.mark.parametrize(('format', 'one_day'), [('cdf_epoch', 86400 * 1000.0),
                                                 ('cdf_epoch16', 86400 + 0j),
                                                 ('cdf_tt2000', 86400 * 10**9)])
def test_parse_time_cdf_units(format, one_day):
    # Check the unit each epoch type counts in: milliseconds, seconds, nanoseconds.
    zero = parse_time(type(one_day)(0), format=format)
    assert (parse_time(one_day, format=format) - zero).sec == pytest.approx(86400.0)


def test_parse_time_cdf_epoch16_keeps_picoseconds():
    # CDF_EPOCH16 exists purely to carry sub-nanosecond times, so the picoseconds
    # in the imaginary part must survive. Casting the complex value to a float
    # would silently truncate this to 2024-01-01T00:00:00.000000000.
    value = complex(CDFepoch.compute([[2024, 1, 1, 0, 0, 0, 123, 456, 789, 12]]))
    assert value == 63871286400 + 123456789012j
    t = parse_time(value, format='cdf_epoch16')
    assert t.utc.isot == '2024-01-01T00:00:00.123456789'
    # astropy's two-double representation resolves a few picoseconds at this
    # epoch, so compare the offset with a tolerance rather than exactly.
    offset = (t - Time(ISO_2024, format='isot', scale='utc')).to(u.ps)
    assert offset.to_value(u.ps) == pytest.approx(123456789012, abs=100)


def test_parse_time_cdf_epoch16_no_complex_warning():
    # Regression test: passing the complex value straight to astropy makes numpy
    # raise a ComplexWarning as it drops the imaginary part, which is easy to miss
    # outside of a test suite that turns warnings into errors.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        warnings.simplefilter("error", np.exceptions.ComplexWarning)
        parse_time(63871286400 + 123456789012j, format='cdf_epoch16')


@pytest.mark.parametrize('nanoseconds', [0, 1, 123, 999999999])
def test_parse_time_cdf_tt2000_is_nanosecond_exact(nanoseconds):
    # CDF_TIME_TT2000 is an int64, which holds more significant digits than a
    # float64: rounding a present-day value through a single float loses ~128 ns.
    t = parse_time(TT2000_2024 + nanoseconds, format='cdf_tt2000')
    assert t.utc.isot == f'2024-01-01T00:00:00.{nanoseconds:09d}'


@pytest.mark.parametrize('nanoseconds', [0, 1, 123, 999999999])
def test_parse_time_cdf_tt2000_before_2000(nanoseconds):
    # TT2000 values before its epoch are negative, where the split into whole
    # seconds plus a remainder has to floor rather than truncate.
    t = parse_time(-315575942816000000 + nanoseconds, format='cdf_tt2000')
    assert t.utc.isot == f'1990-01-01T00:00:00.{nanoseconds:09d}'


@pytest.mark.parametrize('tt2000', [np.iinfo(np.int64).min, np.iinfo(np.int64).max])
def test_parse_time_cdf_tt2000_extremes(tt2000):
    # The ends of the int64 range are what CDF uses for fill and pad values, so
    # they do turn up in real files. They are ~292 years either side of the epoch,
    # and must not wrap around to the wrong side of it.
    t = parse_time(int(tt2000), format='cdf_tt2000')
    assert (t > Time('2000-01-01')) == (tt2000 > 0)
    assert abs((t - Time('2000-01-01')).to_value(u.yr)) == pytest.approx(292, abs=1)


def test_parse_time_cdf_tt2000_leap_second():
    # The last minute of 2016 was 61 seconds long, so UTC ran to 23:59:60. TT2000
    # counts in TT, so unlike a UTC-based epoch it can name that second.
    before, leap, after = (parse_time(v, format='cdf_tt2000') for v in
                           (536500867184000000, 536500868184000000, 536500869184000000))
    assert before.utc.isot == '2016-12-31T23:59:59.000000000'
    assert leap.utc.isot == '2016-12-31T23:59:60.000000000'
    assert after.utc.isot == '2017-01-01T00:00:00.000000000'
    # The three are consecutive seconds, which is what makes the middle one a leap
    # second rather than a duplicate of either neighbour.
    assert (leap - before).sec == pytest.approx(1.0)
    assert (after - leap).sec == pytest.approx(1.0)
    # cdflib agrees on the seconds either side. It writes the leap second itself as
    # "23:60:00", so there is nothing to compare strings against for that one.
    assert CDFepoch.encode(536500867184000000) == '2016-12-31T23:59:59.000000000'
    assert CDFepoch.encode(536500869184000000) == '2017-01-01T00:00:00.000000000'


@pytest.mark.parametrize(('format', 'value', 'scale'), ALL_FORMATS)
def test_parse_time_cdf_agrees_with_cdflib_to_datetime(format, value, scale):
    # cdflib.epochs.CDFepoch.to_datetime is the route most CDF readers take, and
    # is implemented from the epoch definitions rather than via astropy.
    expected = np.datetime64(CDFepoch.to_datetime(value)[0])
    assert parse_time(value, format=format).utc.datetime64 == expected


@pytest.mark.parametrize(('format', 'value', 'scale'), ALL_FORMATS)
def test_parse_time_cdf_default_precision(format, value, scale):
    # The default of 9 digits means nanoseconds are shown rather than dropped.
    assert parse_time(value, format=format).precision == 9


@pytest.mark.parametrize(('format', 'value', 'scale'), ALL_FORMATS)
def test_parse_time_cdf_passes_kwargs_to_time(format, value, scale):
    # Keyword arguments must reach astropy rather than being silently dropped.
    assert parse_time(value, format=format, precision=3).precision == 3
    with pytest.raises(TypeError, match='unexpected keyword argument'):
        parse_time(value, format=format, not_a_time_kwarg=True)


def test_cdf_formats_registered_with_astropy():
    # A user should not have to import cdflib themselves for the formats to work,
    # and once parse_time has run they are usable directly from astropy too.
    parse_time(TT2000_2024, format='cdf_tt2000')
    assert {'cdf_epoch', 'cdf_epoch16', 'cdf_tt2000'} <= set(Time.FORMATS)


@pytest.mark.parametrize(('format', 'value', 'scale'), ALL_FORMATS)
def test_is_time_cdf(format, value, scale):
    # The formats work through the other public entry point as well.
    assert is_time(value, time_format=format)


@pytest.mark.parametrize('format', ['cdf_epoch', 'cdf_epoch16', 'cdf_tt2000'])
def test_parse_time_cdf_without_cdflib(format):
    # Simulate cdflib being absent: a None entry in sys.modules makes the import
    # machinery raise ImportError even though cdflib is installed here.
    with mock.patch.dict(sys.modules, {'cdflib': None, 'cdflib.epochs_astropy': None}):
        with pytest.raises(ImportError, match=f"cdflib must be installed to use the '{format}' time format"):
            parse_time(0, format=format)


def test_parse_time_docstring_mentions_cdf_formats():
    # Users have to know the format names exist before they can pass one.
    assert "'cdf_epoch'" in parse_time.__doc__
    assert "'cdf_epoch16'" in parse_time.__doc__
    assert "'cdf_tt2000'" in parse_time.__doc__
