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
import astropy.units as u

import sunpy.timeseries
from sunpy.data.sample import RHESSI_TIMESERIES

###############################################################################
# Let's first load a TimeSeries from sample data.
# Input data contains 9 columns, which are evenly sampled with a time
# step of 4 seconds.So, the data is suitable for signal analysis.
ts = sunpy.timeseries.TimeSeries(RHESSI_TIMESERIES)

###############################################################################
# We now use SciPy's periodogram to estimate the power spectra of the first 
# column of the Timeseries.First column contains X-Ray emmisions in the range 
# of 3-6 keV.Astropy's Lomb-Scargle Periodograms 
# (http://docs.astropy.org/en/v2.0.5/stats/lombscargle.html) can be used
# instead of SciPy's periodogram.

x_ray = ts.columns[0]
# The suitable value for fs would be 0.25 Hz as the time step is 4 s.
freq, spectra = signal.periodogram(ts.data[x_ray], fs=0.25)

###############################################################################
# Plot the power spectrum

plt.semilogy(freq, spectra)
plt.title('Power Spectrum of {}'.format(x_ray))
plt.ylabel('Power Spectral Density [{:LaTeX}]'.format(ts.units[x_ray] ** 2 / u.Hz))
plt.xlabel('Frequency [Hz]')
plt.show()
