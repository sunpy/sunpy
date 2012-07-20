from datetime import datetime

import numpy as np

from scipy import ndimage

from sunpy.spectra.spectrogram import Spectrogram, LinearTimeSpectrogram


def mk_spec(image):
	return Spectrogram(
		image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
		np.linspace(0, image.shape[0] - 1, image.shape[0]),
		datetime(2010, 10, 10), datetime(2010, 10, 10, 1), 0
	)

def test_subtract_bg():
	bg = np.linspace(0, 200, 200).astype(np.uint16)
	bg.shape = (200, 1)
	bg = bg + np.zeros((200, 3600))

	signal = np.random.rand(200, 1800) * 255
	signal = signal.astype(np.uint16)

	image = bg
	image[:, 1800:] += signal

	spectrogram = mk_spec(image)
	assert np.array_equal(
		spectrogram.subtract_bg()[:, 1800:], signal
	)


def test_slice_time_axis():
	spectrogram = mk_spec(np.zeros((200, 3600)))
	new = spectrogram[:, 59:]
	assert new.t_init == 59
	assert np.array_equal(new.time_axis, np.linspace(0, 3600 - 60, 3600 - 59))
	assert new.start == datetime(2010, 10, 10, 0, 0, 59)


def test_slice_freq_axis():
	spectrogram = mk_spec(np.zeros((200, 3600)))
	new = spectrogram[100:, :]
	assert np.array_equal(new.freq_axis, np.linspace(100, 199, 100))


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

def test_join():
	image = np.random.rand(200, 3600)
	one = LinearTimeSpectrogram(
		image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
		np.linspace(0, image.shape[0] - 1, image.shape[0]),
		datetime(2010, 10, 10), datetime(2010, 10, 10, 0, 30), 0, 0.5,
	)

	image = np.random.rand(200, 3600)
	other = LinearTimeSpectrogram(
		image, np.linspace(0, image.shape[1] - 1, image.shape[1]),
		np.linspace(0, image.shape[0] - 1, image.shape[0]),
		datetime(2010, 10, 10, 0, 29), datetime(2010, 10, 10, 1, 29), 1799, 1,
	)

	z = LinearTimeSpectrogram.join_many([one, other])
	assert z.shape == (200, 3 * 3600 - 2)

	assert np.array_equal(z[:, :3598], one[:, :-2])
	assert np.array_equal(z[:, 3598:], ndimage.zoom(other, (1, 2)))
