"""
===============================================
Parallel Spectral Fitting with Astropy Modeling
===============================================

In the 7.0 release of astropy, a new function was added `~astropy.modeling.fitting.parallel_fit_dask`, alongside multiple performance improvements to fitting astropy models with a non-linear fitter.
These changes mean that it is now practical to fit many independent models to a large multi-dimensional array of data.

In this example we will demonstrate this functionality by fitting many spectra of a raster scan in an observation of the solar chromosphere by the `SPICE <https://spice.ias.u-psud.fr/>`__ instrument on the Solar Orbiter mission.
"""

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ndcube import NDCube

import astropy.units as u
from astropy.io import fits
from astropy.modeling import models as m
from astropy.modeling.fitting import TRFLSQFitter, parallel_fit_dask
from astropy.wcs import WCS

###############################################################################
# For this example, we are going to use a part of an observation from the SPICE instrument
# which is a rastering spectrograph onboard Solar Orbiter. The focus will be on the
# spectral window containing the Nitrogen IV line (76.51 nm) and the Neon VIII line
# (77.04 nm).

res = Fido.search(a.Time("2023-04-15 01:00", "2023-04-15 02:00"), 
                  a.Instrument.spice, a.Level(2), 
                  a.Provider.soar, 
                  a.soar.Product.spice_n_ras)
filename = Fido.fetch(res)[0]

###############################################################################
# We now open the FITS file and access the window via EXTNAME.

hdul = fits.open(filename)
hdu = hdul['N IV 765 ... Ne VIII 770 (Merged)']

###############################################################################
# The next step is to create an `~ndcube.NDCube` object from the data we have opened.
#
# To construct this `~ndcube.NDCube` object we use the header from the
# FITS file to make a WCS object, the data and the unit of the data,
# as well as constructing a mask for all points where the data is NaN.
#
# We then crop down the cube to make it faster to work with.

spice = NDCube(hdu.data, wcs=WCS(hdu), unit=hdu.header["BUNIT"], mask=np.isnan(hdu.data))
# The first dimension is length one so we will drop it
spice = spice[0]
# To ensure this example is quick, we will only do a 50x50 box
spice = spice[:, 100:150, 100:1500]

###############################################################################
# This cube we have constructed has three pixel and four world
# dimensions, the first array dimension is wavelength (71 long)
# followed by the slit dimension and then the rastering dimension. The
# world axes are wavelength, helioprojective latitude and longitude
# and time which increases along the rastering dimension.

print(spice)

###############################################################################
# To aid our analysis we are going to make two rebinned cubes from
# this data, one summed along the wavelength dimension and one of the
# spectra averaged over all spatial pixels.

wl_sum = spice.rebin((spice.data.shape[0], 1, 1), operation=np.sum)[0]
print(wl_sum)

spatial_mean = spice.rebin((1, *spice.data.shape[1:]))[:, 0, 0]
print(spatial_mean)

###############################################################################
# We can now use the `~ndcube.NDCube` built-in plotting routines to show the
# spatially averaged spectra, which uses `~astropy.visualization.wcsaxes.WCSAxes`.

plt.figure()
ax = spatial_mean.plot(axes_units=[u.nm])
ax.coords[0].set_major_formatter("x.xx")

###############################################################################
# Now we can create a model for this spectra.

NIV_wave = 76.51 * u.nm
NeVIII_wave = 77.04 * u.nm

initial_model = (
    m.Const1D(amplitude=0.1 * spice.unit) +
    m.Gaussian1D(amplitude=1 * spice.unit, mean=NIV_wave, stddev=0.05 * u.nm) +
    m.Gaussian1D(amplitude=1 * spice.unit, mean=NeVIII_wave, stddev=0.05 * u.nm)
)
print(initial_model)

###############################################################################
# To improve our initial conditions we now fit the initial model to
# the spatially averaged spectra. To do this we use the
# `~ndcube.NDCube.axis_world_coords` method of `~ndcube.NDCube` which returns
# all, or a subset of, the world coordinates along however many array
# axes they are correlated with. In this case we get the wavelength
# dimension which only returns a single `~astropy.coordinates.SpectralCoord` object
# corresponding to the first array dimension of the cube.

fitter = TRFLSQFitter(calc_uncertainties=True)
average_fit = fitter(
    initial_model,
    spatial_mean.axis_world_coords("em.wl")[0].to(u.nm),
    spatial_mean.data*spatial_mean.unit,
)
print(average_fit)

###############################################################################
# Now we can add to our previous plot the initial model and the model
# fit to the average spectra.

fig = plt.figure()
ax = spatial_mean.plot(label="spatial average")
ax.coords[0].set_format_unit("nm")
ax.coords[0].set_major_formatter("x.xx")
ax.plot(initial_model(spatial_mean.axis_world_coords("em.wl")[0].to(u.AA)), label="initial model")
ax.plot(average_fit(spatial_mean.axis_world_coords("em.wl")[0].to(u.AA)), linestyle="--", label="spatial average fit")
plt.legend()

###############################################################################
# Parallel Fitting
# ----------------

###############################################################################
# Now we have our model to fit to all our spectra we can start working on the parallel fitting.
#
# The function `~astropy.modeling.fitting.parallel_fit_dask` will map
# a model to each element of a cube along one (or more) "fitting
# axes", in this case our fitting axis is our wavelength axis (array
# axis 0). So we want to fit each slice of the data array along the
# 0th axis.
#
# The key arguments to the `~astropy.modeling.fitting.parallel_fit_dask` function are:
#
# * A data array. This can be a numpy array or a dask array, or a `~astropy.nddata.NDData` / `~ndcube.NDCube` object. If it's one of these objects then the data, wcs, mask, data_unit and uncertainty are all extracted from the object and used in place of their respective keyword arguments.
# * A model to fit.
# * A fitter instance.
# * The fitting axis (or axes).
#
# What is returned from `~astropy.modeling.fitting.parallel_fit_dask` is a model with array parameters with the shape of the non-fitting axes of the data (so in this case 100x100 arrays).

###############################################################################
# We can therefore fit all our SPICE cube as follows:
# .. note::
#
#     All examples here are done with ``scheduler="single-threaded"``,
#     this is to allow them to build on our documentation. When running this
#     example yourself, you should remove this line.

spice_model_fit = parallel_fit_dask(
    data=spice,
    model=average_fit,
    fitter=TRFLSQFitter(),
    fitting_axes=0,
    scheduler="single-threaded",
)

###############################################################################
# Given that we are going to want to visualize the output of a few fits,
# we will define a plotting function which will display the shift in the peak
# locations of the two Gaussians. We shall talk more about this later.

def plot_spice_fit(spice_model_fit):
    g1_peak_shift = spice_model_fit.mean_1.quantity.to(u.km/u.s, equivalencies=u.doppler_optical(NIV_wave))
    g2_peak_shift = spice_model_fit.mean_2.quantity.to(u.km/u.s, equivalencies=u.doppler_optical(NeVIII_wave))

    fig, axs = plt.subplots(nrows=3, subplot_kw=dict(projection=wl_sum), figsize=(5,  15))
    fig.suptitle(f"SPICE - {hdu.header["EXTNAME"]} - {hdu.header["DATE-AVG"]}")

    wl_sum.plot(axes=axs[0])
    fig.colorbar(axs[0].get_images()[0], ax=axs[0], extend="both", label=f"{wl_sum.unit:latex}", shrink=0.8)
    axs[0].set_title("Data (summed over wavelength)", pad=40)

    g1_max = np.percentile(np.abs(g1_peak_shift.value), 99)
    mean_1 = axs[1].imshow(g1_peak_shift.value, cmap="coolwarm", vmin=-g1_max, vmax=g1_max)
    fig.colorbar(mean_1, ax=axs[1], extend="both", label=f"Velocity from Doppler shift [{g1_peak_shift.unit:latex}]", shrink=0.8)
    axs[1].set_title(f"N IV ({NIV_wave:latex})", pad=40)

    g2_max = np.percentile(np.abs(g2_peak_shift.value), 98)
    mean_2 = axs[2].imshow(g2_peak_shift.value, cmap="coolwarm", vmin=-g2_max, vmax=g2_max)
    fig.colorbar(mean_2, ax=axs[2], extend="both", label=f"Velocity from Doppler shift [{g2_peak_shift.unit:latex}]", shrink=0.8)
    axs[2].set_title(f"Ne VIII ({NeVIII_wave:latex})", pad=40)

    for ax in axs:
        ax.coords[0].set_ticklabel(exclude_overlapping=True)
        ax.coords[0].set_axislabel("Helioprojective Longitude")
        ax.coords[1].set_axislabel("Helioprojective Latitude")
        ax.coords[2].set_ticklabel(exclude_overlapping=True)

    fig.tight_layout()


plot_spice_fit(spice_model_fit)

###############################################################################
# Oh dear! This clearly didn't work.
#
# To discover why we can use the "diagnostics" functionality of the `~astropy.modeling.fitting.parallel_fit_dask` function.
# This lets each separate process write out logs of errors
# or warnings to a directory of our choice, or run a function (useful for making diagnostic plots).
# In this case we are going to have it write out error logs.

###############################################################################
# First we define the local path we want the logs saved to and ensure the directory
# and the contents of that directory have been removed (to make sure that
# no output from previous runs is present).

diag_path = Path("./diag")
shutil.rmtree(diag_path, ignore_errors=True)

###############################################################################
# We pass the ``diagnostics="error"`` argument to enable logging of error messages
# and the ``diagnostics_path=`` argument to specify where to save the logs.

spice_model_fit = parallel_fit_dask(
    data=spice,
    model=average_fit,
    fitter=TRFLSQFitter(),
    fitting_axes=0,
    diagnostics="error",
    diagnostics_path=diag_path,
    scheduler="single-threaded",
)

###############################################################################
# We can now find all the folders in the diagnostics path:

diag_folders = list(diag_path.glob("*"))

###############################################################################
# And read the contents of each log into a list.

errors = []
for diag in diag_folders:
    if (path := (diag/"error.log")).exists():
        content = open(path).read()
        errors.append(content)

###############################################################################
# We can now print out the first error.

print(f"{len(errors)} errors occurred")
print("First error is:")
print(errors[0])

###############################################################################
# The reason for the failure of the fitting was the presence of NaN values in the data array.
# When calling the `~astropy.modeling.fitting.TRFLSQFitter` (and many others) it's possible to set the ``filter_non_finite=True`` keyword argument.
# To do this with the ``parallel_fit_dask`` function we pass a dictionary of keyword arguments to the fitter as the ``fitter_kwargs`` argument:

spice_model_fit = parallel_fit_dask(
    data=spice,
    model=average_fit,
    fitter=TRFLSQFitter(),
    fitting_axes=0,
    fitter_kwargs={"filter_non_finite": True}, # Filter out non-finite values,
    scheduler="single-threaded",
)

plot_spice_fit(spice_model_fit)

###############################################################################
# That's better!

###############################################################################
# Working with the fit
# --------------------

###############################################################################
# The return value of the `~astropy.modeling.fitting.parallel_fit_dask` function
# is an Astropy model instance with the parameters set based on the result of the fit.
# This is the same as the return value of the fitter called in serial, so for more
# information about how to work with the results of the fit, you can read the
# Astropy documentation for serial fitting, such as this page on
# `Fitting Models to Data <https://docs.astropy.org/en/stable/modeling/fitting.html>`__.
#
# We shall quickly cover some key points.

###############################################################################
# Our input model is a ``CompoundModel`` which is a model that combines many models
# together via various operators, in our case the ``+`` operator. Each individual model
# can be accessed by using slicing notation, if we print out our model we can see that
# we have three models (0, 1 and 2) all added together.

print(spice_model_fit)

###############################################################################
# To access the first Gaussian model we can do this:

print(spice_model_fit[1])

###############################################################################
# If we want to access the parameters of this model we can do it in two different ways:

print(spice_model_fit[1].mean)

###############################################################################
# or
print(spice_model_fit.mean_1)

###############################################################################
# In our plotting helper above we access the mean parameters of both the Gaussian
# fits, let's take a closer look at that.
# The parameters on a model are `~astropy.modeling.Parameter` classes,
# but they can be converted to `~astropy.units.Quantity` objects by accessing
# their ``.quantity`` property:

print(spice_model_fit.mean_1.quantity)

###############################################################################
# A `~astropy.units.Quantity` object can be converted to other units:

print(spice_model_fit.mean_1.quantity.to(u.AA))

###############################################################################
# Using :ref:`unit_equivalencies` you can do unit conversions which require an
# assumption or some extra calculation. Some of the built-in equivalencies in
# Astropy are for doppler shifts, we can use the `~astropy.units.doppler_optical`
# equivalency to convert to velocity.

spice_model_fit.mean_1.quantity.to(u.km/u.s, equivalencies=u.doppler_optical(NIV_wave))

###############################################################################
# One other thing we may want to do is to evaluate the fitted model for all pixels,
# for example, to plot them or otherwise inspect a single fit.
# This can be done by passing in a wavelength array which is `broadcastable <https://numpy.org/doc/stable/user/basics.broadcasting.html>`__
# to the shape of the non-fitting axes. We can do this by once again using the
# `~ndcube.NDCube.axis_world_coords` method of `~ndcube.NDCube`.

wavelength = spatial_mean.axis_world_coords("em.wl")[0]

print(wavelength.shape)

###############################################################################
# And then to make it the correct shape we can add two dummy dimensions to the end:

wavelength = wavelength[:, None, None]

print(wavelength.shape)

###############################################################################
# We can now evaluate the model for all spatial points with this input:

all_fits = spice_model_fit(wavelength)

###############################################################################
# We can now make a plot of all the fits.

fig = plt.figure(figsize=(11, 5))

ax = spatial_mean.plot(axes_units=[u.nm], label="Average spectra", zorder=99)
ax.coords["wavelength"].set_major_formatter("x.xx")

ax.plot(average_fit(spatial_mean.axis_world_coords("em.wl")[0]),
        linestyle="--",
        label="Spatial Average Model",
        zorder=99)

# Iterate over each fit and plot it.
for fit_arr in all_fits.reshape((spatial_mean.data.shape[0], -1)).T:
    line, = ax.plot(fit_arr, alpha=0.3, color="C3", linewidth=0.05)

# Use the last line plotted as the anchor for the legend
line.set_label("Fitted Pixels")
leg = plt.legend(loc="upper left")
# Make the legend not semi-transparent
for lh in leg.legend_handles:
    lh.set_alpha(1)

plt.show()
