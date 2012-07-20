from datetime import datetime

import numpy as np

from sunpy.spectra.spectrogram import Spectrogram

def test_subtract_bg():
	bg = np.linspace(0, 200, 200).astype(np.uint16)
	bg.shape = (200, 1)
	bg = bg + np.zeros((200, 3600))

	signal = np.random.rand(200, 1800) * 255
	signal = signal.astype(np.uint16)

	image = bg
	image[:, 1800:] += signal

	spectrogram = Spectrogram(
		image, np.linspace(0, 3599, 3600), np.linspace(0, 199, 200),
		datetime(2010, 10, 10), datetime(2010, 10, 10, 1), 5
	)
	assert np.array_equal(
		spectrogram.subtract_bg()[:, 1800:], signal
	)

