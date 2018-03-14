"""
===============================
Power Spectrum of a TimeSeries
===============================

An example showing how to estimate the power spectrum of a TimeSeries.
"""

##############################################################################
# Start by importing the necessary modules.

import matplotlib.pyplot as plt
import pandas as pd 
from scipy import signal

import sunpy.timeseries
from sunpy.data.sample import GOES_XRS_TIMESERIES

###############################################################################
# Let's first load a TimeSeries from sample data.

ts = sunpy.timeseries.TimeSeries(GOES_XRS_TIMESERIES)
# Data is taken from GOES X-Ray Sensor(XRS), which observes full-disk-integrated
# solar flux in two broadband channels: 1--8 angstrom (long); and
# 0.5--4 angstrom (short)

###############################################################################
# We now use scipy's periodogram to estimate the power spectra of the first 
# column of the Timeseries

xray_short = ts.columns[0]
freq, spectra = signal.periodogram(ts.data[xray_short], fs=1)

###############################################################################
# Plot the power spectrum

plt.semilogy(freq, spectra)
plt.title('Power Spectrum of {}'.format(xray_short))
plt.ylabel('Power Spectral Density(({})**2/Hz)'.format(ts.units[xray_short]))
plt.xlabel('Frequency(Hz)')
plt.show()





