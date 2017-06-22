"""
=============
Map Histogram
=============

How to plot and fit the histogram of the data of a map.
"""
import numpy as np
import matplotlib.pyplot as plt
from sunpy.data.sample import AIA_171_IMAGE
from sunpy.map import Map
from astropy.modeling import models, fitting

###############################################################################
# We first create the Map using the sample data.
aia = Map(AIA_171_IMAGE)

###############################################################################
# We will fit the log of the data so first lets take the log and then get the
# histogram of that data. Since we cannot log negative values we remove them.
data = np.log10(aia.data[aia.data > 0])
hist, bins = np.histogram(data, bins=np.linspace(0,10,30))
width = 0.7 * (bins[1] - bins[0])
x = (bins[:-1] + bins[1:]) / 2

###############################################################################
# Next we will fit the histogram using astropy models which provide a simple
# interface to fitting functions.
g_init = models.Gaussian1D(amplitude=hist.max(), mean=data.mean(), stddev=data.std())
fit = fitting.LevMarLSQFitter()
fitted_model = fit(g_init, x, hist)

###############################################################################
# Finally let's plot the histogram and the fit.
plt.figure()
fit_params_str = "height={0:.2f}\n".format(fitted_model.amplitude.value)
fit_params_str += "mean={1:.2f}\n".format(fitted_model.mean.value)
fit_params_str += "sigma={2:.2f}".format(fitted_model.stddev.value)
plt.bar(x, hist, align='center', width=width, label='Histogram')
plt.plot(x, fitted_model(x), 'r-', label="Gaussian Fit")
plt.xlim(0, 10)
plt.xlabel('Intensity of LOG10(data)')
plt.text(5, 150000, fit_params_str)
plt.legend()
plt.savefig('gallery_histogram.png')
