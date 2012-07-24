# -*- coding: utf-8 -*-
# Author: Florian Mayer <florian.mayer@bitsrc.org>

from __future__ import absolute_import

from datetime import datetime

import numpy as np

from sunpy.data.sample import CALLISTO_IMAGE
from sunpy.spectra.sources.callisto import CallistoSpectrogram


def test_read():
	ca = CallistoSpectrogram.read(CALLISTO_IMAGE)
	assert ca.start == datetime(2011, 9, 22, 10, 30, 0, 51000)
	assert (
		ca.t_init ==
		(datetime(2011, 9, 22, 10, 30) - datetime(2011, 9, 22)).seconds
	)
	assert ca.shape == (200, 3600)
	assert ca.t_delt == 0.25
	# Test linearity of time axis.
	assert np.array_equal(
        ca.time_axis, np.linspace(0, 0.25 * (ca.shape[1] - 1), ca.shape[1])
    )
