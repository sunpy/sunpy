"""
======================================================
Smoothing of TimeSeries Data Using Convolution Filters
======================================================

Sometimes, when working with scientific data, we have noisy data from which
we need to extract low-frequency components from. This example illustrates
creating a TimeSeries data, adding noise to it and then smoothing it using
convolution filter kernels from `~astropy.convolution` and 
`~astropy.convolution.convolve` function.
"""

##############################################################################
# Start by importing the necessary modules.

import datetime

import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame

from astropy.convolution import convolve, Box1DKernel

from sunpy.timeseries import TimeSeries

###############################################################################
# We will first create a sine wave with time range spanning over 1 day and then
# convert it into a `TimeSeries` object by passing as pandas DataFrame.

base = datetime.datetime.today()
times = [base - datetime.timedelta(minutes=x) for x in range(0, 24 * 60)]
intensity = np.sin(np.arange(0, 12 * np.pi, step=(12 * np.pi) / (24 * 60)))
original_df = DataFrame(intensity, index=times, columns=['intensity'])
# Convert into TimeSeries
original_ts = TimeSeries(original_df)
original_ts.peek(title="Original TimeSeries")

###############################################################################
# Now we will create a noise that is to be added to original wave and convert 
# resulting signal into a `TimeSeries` object.

noise = (np.random.rand(24 * 60) - 0.5)
noisy_signal = intensity + noise
noisy_signal_df = DataFrame(noisy_signal, index=times, columns=['intensity'])
noisy_signal_ts = TimeSeries(noisy_signal_df)
noisy_signal_ts.peek(title="Noisy TimeSeries")

###############################################################################
# We can now apply the convolution filter on noisy signal to get smooth data.

data = noisy_signal_ts.data['intensity'].values
# Apply convolution filter
convolved_data = convolve(data, kernel=Box1DKernel(20))
convolved_data_df = DataFrame(convolved_data, index=times, columns=['intensity'])
convolved_data_ts = TimeSeries(convolved_data_df)
convolved_data_ts.peek(title="Smoothed TimeSeries")
