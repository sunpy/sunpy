# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

from datetime import datetime
import pytest

import numpy as np

from numpy.testing import assert_array_almost_equal

from sunpy.spectra.spectrogram import(
    Spectrogram, LinearTimeSpectrogram, _LinearView
)


def is_linear(arr):
    return np.array_equal(arr, np.linspace(arr[0], arr[-1], len(arr)))


def dict_eq(one, other):
    ks = set(one.keys())
    if ks != set(other.keys()):
        return False
    for key in ks:
        if isinstance(one[key], np.ndarray):
            if not np.array_equal(one[key], other[key]):
                return False
        else:
            if one[key] != other[key]:
                return False
    return True


def mk_spec(image):
    return Spectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10), datetime(2010, 10, 10, 1), 0
    )


def test_subtract_bg():
    # The idea is to generate background and add a random signal, perform
    # background subtraction and see if the signal comes out again.
    bg = np.linspace(0, 200, 200).astype(np.uint16)
    bg.shape = (200, 1)
    bg = bg + np.zeros((200, 3600))

    signal = np.random.rand(200, 1800) * 255
    signal = signal.astype(np.uint16)

    image = bg
    image[:, 1800:] += signal

    spectrogram = mk_spec(image)
    sbg = spectrogram.subtract_bg()
    assert np.array_equal(
        spectrogram.subtract_bg()[:, 1800:].data, signal
    )

    assert dict_eq(spectrogram._get_params(), sbg._get_params())


def test_auto_const_bg():
    # The idea is to generate background and add a random signal, perform
    # background subtraction and see if the signal comes out again.
    x = np.linspace(0, 200, 200).astype(np.uint16)
    bg = x.reshape(200, 1)
    bg = bg + np.zeros((200, 3600))

    signal = np.random.rand(200, 1800) * 255
    signal = signal.astype(np.uint16)

    image = bg
    image[:, 1800:] += signal

    spectrogram = mk_spec(image)
    sbg = spectrogram.auto_const_bg()
    assert np.array_equal(sbg, x.reshape(200, 1))


def test_randomized_auto_const_bg():
    # The idea is to generate background and add a random signal, perform
    # background subtraction and see if the signal comes out again.
    # As this is a Monte-Carlo probabilistic algorithm this test might
    # fail occasionally.
    x = np.linspace(0, 200, 200).astype(np.uint16)
    bg = x.reshape(200, 1)
    bg = bg + np.zeros((200, 3600))

    signal = np.random.rand(200, 1800) * 255
    signal = signal.astype(np.uint16)

    image = bg
    image[:, 1800:] += signal

    spectrogram = mk_spec(image)
    sbg = spectrogram.randomized_auto_const_bg(1500)
    assert np.array_equal(sbg, x.reshape(200, 1))


def test_slice_time_axis():
    rnd = np.random.rand(200, 3600)
    spectrogram = mk_spec(rnd)
    new = spectrogram[:, 59:3599]
    assert new.shape == (200, 3600 - 59 - 1)
    assert new.t_init == 59
    assert np.array_equal(new.time_axis,
        np.linspace(0, 3600 - 60 - 1, 3600 - 59 - 1)
    )
    assert new.start == datetime(2010, 10, 10, 0, 0, 59)
    assert np.array_equal(new.data, rnd[:, 59:3599])


def test_slice_freq_axis():
    rnd = np.random.rand(200, 3600)
    spectrogram = mk_spec(rnd)
    new = spectrogram[100:150, :]
    assert new.shape == (50, 3600)
    assert np.array_equal(new.freq_axis, np.linspace(100, 149, 50))
    assert np.array_equal(new.data, rnd[100:150, :])


def test_slice_both_axis():
    rnd = np.random.rand(200, 3600)
    spectrogram = mk_spec(rnd)
    new = spectrogram[100:, 59:]
    assert new.shape == (100, 3600 - 59)
    assert new.t_init == 59
    assert np.array_equal(new.time_axis, np.linspace(0, 3600 - 60, 3600 - 59))
    assert new.start == datetime(2010, 10, 10, 0, 0, 59)
    assert np.array_equal(new.freq_axis, np.linspace(100, 199, 100))
    assert np.array_equal(new.data, rnd[100:, 59:])


def test_time_to_x():
    image = np.zeros((200, 3600))
    spectrogram = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10), datetime(2010, 10, 10, 1), 0, 1
    )
    ret = spectrogram.time_to_x(datetime(2010, 10, 10, 0, 0, 59))
    assert isinstance(ret, int)
    assert ret == 59


def test_time_to_x_nonlinear():
    image = np.zeros((200, 3600))
    spectrogram = Spectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10), datetime(2010, 10, 10, 1)
    )
    ret = spectrogram.time_to_x(datetime(2010, 10, 10, 0, 0, 59))
    assert isinstance(ret, int)
    assert ret == 59


def test_join():
    image = np.random.rand(200, 3600)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10), datetime(2010, 10, 10, 0, 30), 0, 0.5,
    )

    image = np.random.rand(200, 3600)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10, 0, 29),
        datetime(2010, 10, 10, 1, 29), 1799, 1,
    )

    z = LinearTimeSpectrogram.join_many(
        [one, other], nonlinear=False, maxgap=0
    )
    # The - 1 is because resampling other produces an image of size
    # 2 * 3600 - 1
    # The - 2 is because there is one second overlap.
    assert z.shape == (200, 3 * 3600 - 2 - 1)

    assert np.array_equal(z.data[:, :3598], one.data[:, :-2])
    # assert np.array_equal(z[:, 3598:], ndimage.zoom(other, (1, 2)))
    assert z.start == one.start
    assert z.end == other.end
    assert is_linear(z.time_axis)
    assert isinstance(z, LinearTimeSpectrogram)


def test_join_dtype():
    image = np.random.rand(200, 3600).astype(np.uint8)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10), datetime(2010, 10, 10, 0, 30), 0, 0.5,
    )

    image = np.random.rand(200, 3600).astype(np.uint8)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10, 0, 29),
        datetime(2010, 10, 10, 1, 29), 1799, 1,
    )

    z = LinearTimeSpectrogram.join_many(
        [one, other], nonlinear=False, maxgap=0
    )
    assert z.dtype == np.dtype('uint8')


def test_join_different_dtype():
    image = np.random.rand(200, 3600).astype(np.uint16)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10), datetime(2010, 10, 10, 0, 30), 0, 0.5,
    )

    image = np.random.rand(200, 3600).astype(np.uint8)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10, 0, 29),
        datetime(2010, 10, 10, 1, 29), 1799, 1,
    )

    z = LinearTimeSpectrogram.join_many(
        [one, other], nonlinear=False, maxgap=0
    )
    assert z.dtype == np.dtype('uint16')


def test_join_midnight():
    image = np.random.rand(200, 3600)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10, 23, 30),
        datetime(2010, 10, 10, 23, 59, 59), 84600, 0.5,
    )

    image = np.random.rand(200, 3600)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 11, 0, 0), datetime(2010, 10, 11, 1), 0, 1,
    )

    z = LinearTimeSpectrogram.join_many(
        [one, other], nonlinear=False, maxgap=0
    )
    # The - 1 is because resampling other produces an image of size
    # 2 * 3600 - 1
    assert z.shape == (200, 3 * 3600 - 1)

    assert np.array_equal(z.data[:, :3600], one.data)
    assert is_linear(z.time_axis)
    assert isinstance(z, LinearTimeSpectrogram)


def test_join_month():
    image = np.random.rand(200, 3600)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2012, 7, 31, 23, 30),
        datetime(2012, 7, 31, 23, 59, 59), 84600, 0.5,
    )

    image = np.random.rand(200, 3600)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2012, 8, 1), datetime(2012, 8, 1, 1), 0, 1,
    )

    z = LinearTimeSpectrogram.join_many(
        [one, other], nonlinear=False, maxgap=0
    )
    # The - 1 is because resampling other produces an image of size
    # 2 * 3600 - 1
    assert z.shape == (200, 3 * 3600 - 1)

    assert np.array_equal(z.data[:, :3600], one.data)
    assert is_linear(z.time_axis)
    assert isinstance(z, LinearTimeSpectrogram)


def test_join_year():
    image = np.random.rand(200, 3600)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2012, 12, 31, 23, 30),
        datetime(2013, 1, 1, 0, 0, 0), 84600, 0.5,
    )

    image = np.random.rand(200, 3600)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2013, 1, 1), datetime(2013, 1, 1, 1), 0, 1,
    )

    z = LinearTimeSpectrogram.join_many(
        [one, other], nonlinear=False, maxgap=0
    )
    # The - 1 is because resampling other produces an image of size
    # 2 * 3600 - 1
    assert z.shape == (200, 3 * 3600 - 1)

    assert np.array_equal(z.data[:, :3600], one.data)
    assert is_linear(z.time_axis)
    assert isinstance(z, LinearTimeSpectrogram)


def test_join_over_midnight():
    image = np.random.rand(200, 3600)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10, 23, 45),
        datetime(2010, 10, 11, 0, 15,), 85500, 0.5,
    )
    image = np.random.rand(200, 3600)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 11, 0, 15), datetime(2010, 10, 11, 1, 15), 900, 1,
    )

    z = LinearTimeSpectrogram.join_many(
        [one, other], nonlinear=False, maxgap=0
    )
    # FIXME: not used?!
    oz = other.resample_time(0.5)

    # The - 1 is because resampling other produces an image of size
    # 2 * 3600 - 1
    assert z.shape == (200, 3 * 3600 - 1)

    assert np.array_equal(z.data[:, :3600], one.data)
    assert np.array_equal(z.time_axis[:3600], one.time_axis)
    assert is_linear(z.time_axis)
    assert isinstance(z, LinearTimeSpectrogram)


def test_join_gap():
    image = np.random.rand(200, 3600)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10, 23, 45),
        datetime(2010, 10, 11, 0, 15,), 85500, 0.5,
    )

    image = np.random.rand(200, 3600)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 11, 0, 15, 1),
        datetime(2010, 10, 11, 1, 15), 901, 1,
    )
    with pytest.raises(ValueError) as excinfo:
        LinearTimeSpectrogram.join_many(
            [one, other], nonlinear=False, maxgap=0
        )

    assert excinfo.value.message == "Too large gap."


def test_join_with_gap():
    image = np.random.rand(200, 3600)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10, 23, 45),
        datetime(2010, 10, 11, 0, 15,), 85500, 0.5,
    )

    image = np.random.rand(200, 3600)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 11, 0, 15), datetime(2010, 10, 11, 1, 15), 901, 1,
    )

    z = LinearTimeSpectrogram.join_many(
        [one, other], nonlinear=False, maxgap=1, fill=0
    )

    # The - 1 is because resampling other produces an image of size
    # 2 * 3600 - 1
    # The + 2 is because there is one second without data inserted.
    assert z.shape == (200, 3 * 3600 + 2 - 1)

    assert np.array_equal(z.data[:, :3600], one.data)
    # Second data to unpack masked array
    assert (z.data.data[:, 3600:3602] == 0).all()
    assert is_linear(z.time_axis)
    assert isinstance(z, LinearTimeSpectrogram)


def test_join_with_gap_fill():
    image = np.random.rand(200, 3600)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10, 23, 45),
        datetime(2010, 10, 11, 0, 15,), 85500, 0.5,
    )

    image = np.random.rand(200, 3600)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 11, 0, 15), datetime(2010, 10, 11, 1, 15), 901, 1,
    )

    z = LinearTimeSpectrogram.join_many(
        [one, other], nonlinear=False, maxgap=2, fill=np.NaN
    )
    # The - 1 is because resampling other produces an image of size
    # 2 * 3600 - 1
    # The + 2 is because there is one second without data inserted.
    assert z.shape == (200, 3 * 3600 + 2 - 1)

    assert np.array_equal(z.data[:, :3600], one.data)

    print type(z.data)

    # Second data to unpack masked array
    assert np.isnan(z.data.data[:, 3600:3602]).all()
    assert is_linear(z.time_axis)
    assert isinstance(z, LinearTimeSpectrogram)


def test_join_nonlinear():
    image = np.random.rand(200, 3600)
    one = LinearTimeSpectrogram(
        image, np.linspace(0, 0.5 * (image.shape[1] - 1), image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 10, 23, 45),
        datetime(2010, 10, 11, 0, 15,), 85500, 0.5,
    )

    image = np.random.rand(200, 3600)
    other = LinearTimeSpectrogram(
        image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 10, 11, 0, 15),
        datetime(2010, 10, 11, 1, 15), 901, 1,
    )

    oz = other.resample_time(0.5)

    z = LinearTimeSpectrogram.join_many(
        [one, other], nonlinear=True, maxgap=2
    )

    # The - 1 is because resampling other produces an image of size
    # 2 * 3600 - 1
    assert z.shape == (200, 3 * 3600 - 1)

    assert np.array_equal(z.data[:, :3600], one.data)
    assert np.array_equal(z.time_axis[:3600], one.time_axis)
    assert np.array_equal(z.time_axis[3600:], oz.time_axis + 1801)
    assert isinstance(z, Spectrogram)


def test_auto_t_init():
    image = np.random.rand(200, 3600)
    assert Spectrogram(image,
        np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30)
    ).t_init == 900


def test_rescale():
    image = np.random.rand(200, 3600) * 43
    spec = Spectrogram(image,
        np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30)
    )

    nspec = spec.rescale()

    assert dict_eq(spec._get_params(), nspec._get_params())
    assert_array_almost_equal(nspec.data.max(), 1)
    assert nspec.data.min() == 0


def test_rescale_error():
    image = np.zeros((200, 3600))
    spec = Spectrogram(image,
        np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30)
    )

    with pytest.raises(ValueError) as excinfo:
        spec.rescale(0, 1)
    assert (
        excinfo.value.message ==
        "Spectrogram needs to contain distinct values."
    )


def test_rescale_error2():
    image = np.random.rand(200, 3600) * 43
    spec = Spectrogram(image,
        np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.linspace(0, image.shape[0] - 1, image.shape[0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30)
    )

    with pytest.raises(ValueError) as excinfo:
        spec.rescale(1, 1)
    assert excinfo.value.message == "Maximum and minimum must be different."


def test_resample():
    image = np.array([[0, 1, 2], [0, 1, 2]])
    spec = LinearTimeSpectrogram(
        image, np.array([0, 1, 2]), np.array([0]),
        datetime(2012, 1, 1), datetime(2012, 1, 1, 0, 0, 3),
        0, 1
    )
    r = spec.resample_time(0.5)
    assert r.shape[1] == 5
    assert np.array_equal(r.time_axis, np.linspace(0, 2, 5))


def test_upsample():
    image = np.array([[0, 1, 2, 3], [0, 1, 2, 3]])
    spec = LinearTimeSpectrogram(
        image, np.array([0, 1, 2]), np.array([0]),
        datetime(2012, 1, 1), datetime(2012, 1, 1, 0, 0, 3),
        0, 1
    )
    r = spec.resample_time(2)
    assert r.shape[1] == 2


def test_combine_freqs():
    image = np.random.rand(5, 3600)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.array([8, 6, 4, 2, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        0.25
    )
    image = np.random.rand(5, 3600)
    spec2 = LinearTimeSpectrogram(image,
        np.linspace(0, image.shape[1] - 1, image.shape[1]),
        np.array([9, 7, 5, 3, 1]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        0.25
    )
    comb = LinearTimeSpectrogram.combine_frequencies([spec, spec2])
    stuff = [spec, spec2]

    # print comb

    for freq in xrange(10):
        assert np.array_equal(
            comb[9 - freq, :], stuff[freq % 2][4 - freq // 2, :]
        )


def test_join_diff_freq():
    image = np.random.rand(5, 3600)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 0.25 * (image.shape[1] - 1), image.shape[1]),
        np.array([8, 6, 4, 2, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        0.25
    )
    image = np.random.rand(5, 3600)
    spec2 = LinearTimeSpectrogram(image,
        np.linspace(0, 0.25 * (image.shape[1] - 1), image.shape[1]),
        np.array([9, 7, 5, 3, 1]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        1800,
        0.25
    )

    with pytest.raises(ValueError) as excinfo:
        LinearTimeSpectrogram.join_many([spec, spec2])
    assert excinfo.value.message == "Frequency channels do not match."


def test_intersect_time():
    image = np.random.rand(5, 3600)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 0.25 * (image.shape[1] - 1), image.shape[1]),
        np.array([8, 6, 4, 2, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        0.25
    )
    image = np.random.rand(5, 3600)
    spec2 = LinearTimeSpectrogram(image,
        np.linspace(0, 0.25 * (image.shape[1] - 1), image.shape[1]),
        np.array([9, 7, 5, 3, 1]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        901,
        0.25
    )

    one, other = LinearTimeSpectrogram.intersect_time(
        [spec, spec2]
    )

    assert one.shape[1] == other.shape[1]
    assert one.shape[1] == 3596
    assert np.array_equal(one.data, spec.data[:, 4:])
    assert np.array_equal(other.data, spec2.data[:, :-4])

    assert np.array_equal(one.time_axis, other.time_axis)
    assert one.t_init == other.t_init
    assert is_linear(one.time_axis)
    assert is_linear(other.time_axis)

def test_check_linearity():
    image = np.random.rand(5, 3600)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 0.25 * (image.shape[1] - 1), image.shape[1]),
        np.array([8, 6, 4, 2, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        0.25
    )
    assert spec.check_linearity()
    spec.time_axis[1] += 0.1
    assert not spec.check_linearity()
    assert spec.check_linearity(0.1)
    spec.time_axis[1] -= 0.1
    # The average stays (almost) the same because there are 3600 items.
    spec.time_axis[1] += 0.2 * 0.25
    assert spec.check_linearity(None, 0.2)


def test_flatten():
    flat = np.arange(5 * 3600)
    image = flat.reshape(5, 3600)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 0.25 * (image.shape[1] - 1), image.shape[1]),
        np.array([8, 6, 4, 2, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        0.25
    )
    assert np.array_equal(flat, spec.data.flatten())


def test_in_interval():
    image = np.random.rand(5, 900)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 1 * (image.shape[1] - 1), image.shape[1]),
        np.array([8, 6, 4, 2, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        1
    )

    assert np.array_equal(spec.in_interval("00:15", "00:30").data, spec.data)


def test_in_interval2():
    image = np.random.rand(5, 900)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 1 * (image.shape[1] - 1), image.shape[1]),
        np.array([8, 6, 4, 2, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        1
    )

    assert np.array_equal(
        spec.in_interval("2010-01-01T00:15:00", "00:30").data, spec.data
    )


def test_linearize():
    image = np.random.rand(5, 900)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 1 * (image.shape[1] - 1), image.shape[1]),
        np.array([20, 10, 5, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        1
    )
    # 0   1   2   3   4   5  6  7  8
    # -------- ----------- ----- ---
    # 20 17.5 15 12.5 10 7.5 5 2.5 0

    linear = spec.linearize_freqs()
    assert ((linear.freq_axis[:-1] - linear.freq_axis[1:]) == 2.5).all()

    assert (linear[0] == image[0, :]).all()
    assert (linear[1] == image[0, :]).all()
    assert (linear[2] == image[0, :]).all()
    assert (linear[3] == image[1, :]).all()
    assert (linear[4] == image[1, :]).all()
    assert (linear[5] == image[1, :]).all()
    assert (linear[6] == image[2, :]).all()
    assert (linear[7] == image[2, :]).all()
    assert (linear[8] == image[3, :]).all()


def test_linear_view():
    image = np.random.rand(5, 900)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 1 * (image.shape[1] - 1), image.shape[1]),
        np.array([20, 10, 5, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        1
    )

    linear = _LinearView(spec)
    # assert ((linear.freq_axis[:-1] - linear.freq_axis[1:]) == 2.5).all()

    assert (linear[0] == image[0, :]).all()
    assert (linear[1] == image[0, :]).all()
    assert (linear[2] == image[0, :]).all()
    assert (linear[3] == image[1, :]).all()
    assert (linear[4] == image[1, :]).all()
    assert (linear[5] == image[1, :]).all()
    assert (linear[6] == image[2, :]).all()
    assert (linear[7] == image[2, :]).all()
    assert (linear[8] == image[3, :]).all()


def test_linear_view_indexerror():
    image = np.random.rand(5, 900)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 1 * (image.shape[1] - 1), image.shape[1]),
        np.array([20, 10, 5, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        1
    )

    linear = _LinearView(spec)
    # assert ((linear.freq_axis[:-1] - linear.freq_axis[1:]) == 2.5).all()
    with pytest.raises(IndexError):
        linear[9]


def test_linear_view_negative():
    image = np.random.rand(5, 900)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 1 * (image.shape[1] - 1), image.shape[1]),
        np.array([20, 10, 5, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        1
    )

    linear = _LinearView(spec)
    # assert ((linear.freq_axis[:-1] - linear.freq_axis[1:]) == 2.5).all()
    assert (linear[8] == image[3, :]).all()
    assert (linear[-1] == image[3, :]).all()


def test_linear_view_freqs():
    image = np.random.rand(5, 900)
    spec = LinearTimeSpectrogram(image,
        np.linspace(0, 1 * (image.shape[1] - 1), image.shape[1]),
        np.array([20, 10, 5, 0]),
        datetime(2010, 1, 1, 0, 15),
        datetime(2010, 1, 1, 0, 30),
        900,
        1
    )

    linear = _LinearView(spec)
    # assert ((linear.freq_axis[:-1] - linear.freq_axis[1:]) == 2.5).all()

    assert linear.get_freq(0) == 20
    assert linear.get_freq(1) == 20
    assert linear.get_freq(2) == 20
    assert linear.get_freq(3) == 10
    assert linear.get_freq(4) == 10
    assert linear.get_freq(5) == 10
    assert linear.get_freq(6) == 5
    assert linear.get_freq(7) == 5
    assert linear.get_freq(8) == 0
