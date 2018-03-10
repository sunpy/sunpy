"""
===============================
POWER SPECTRUM OF A TIMESERIES
===============================

An example showing how to estimate power spectrum of a TimeSeries
"""

##############################################################################
# Start by importing the necessary modules.

import matplotlib.pyplot as plt
import pandas as pd 
from scipy import signal

import sunpy.timeseries
from sunpy.data.sample import GOES_XRS_TIMESERIES

###############################################################################
# Let's first load a TimeSeries from sample data

ts = sunpy.timeseries.TimeSeries(GOES_XRS_TIMESERIES)

###############################################################################
# We now use scipy's periodogram to estimate the power spectra of the first 
# column of the Timeseries

cols=ts.columns
freq, spectra = signal.periodogram(ts.data[cols[0]], fs=1)

###############################################################################
# Plot the power spectrum

plt.semilogy(freq, spectra)
plt.title('Power Spectrum of {}'.format(cols[0]))
plt.ylabel('Power Spectral Density')
plt.xlabel('Frequency(Hz)')
plt.show()





