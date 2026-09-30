"""
===============================================
Parallel Spectral Fitting with Astropy Modeling
===============================================

This example demonstrates how to fit a model independently to every spectrum in a spectral cube using `~astropy.modeling.fitting.parallel_fit_dask`.
The technique here can be applied to any spectral cube, but here we use a specific example dataset of a raster scan from the `SPICE <https://spice.ias.u-psud.fr/>`__ instrument onboard Solar Orbiter, focusing on the O VI 1032 Å line.

The SPICE instrument on Solar Orbiter is an EUV imaging spectrograph that builds up rasters of the solar transition region and low corona by stepping its slit across the Sun.
Each raster pixel contains a full spectrum, so a single observation can produce tens of thousands of spectra that we'd like to fit independently, for example, to derive a Doppler velocity map from a line centroid.

This example makes use of the `~astropy.modeling.fitting.parallel_fit_dask` helper and multiple performance improvements added in astropy 7.0 for fitting astropy models with non-linear fitters.
Together these changes mean it is now practical to fit many independent models to a large multi-dimensional array of data.
For more information on using this functionality, refer to the :ref:`astropy:parallel-fitting` page in the astropy documentation.
This example demonstrates the *fitting workflow* and is **not** a science-grade Doppler analysis.
To keep the focus on the fitting, we use the SPICE Level 2 data as-is and skip several corrections that a real science analysis would need to apply to the SPICE data.
"""
# sphinx_gallery_tags = ["Solar Orbiter", "SPICE", "Spectral Fitting", "Visualization"]
# sphinx_gallery_thumbnail_number = -2
# sphinx_gallery_multi_image = "single"

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import CenteredNorm, LogNorm
from ndcube import NDCube

import astropy.units as u
from astropy.io import fits
from astropy.modeling import models as m
from astropy.modeling.fitting import TRFLSQFitter, parallel_fit_dask
from astropy.wcs import WCS

from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# Downloading the Data
# --------------------
#
# For this example, we are going to use a part of an observation from
# the SPICE instrument which is a rastering spectrograph onboard Solar
# Orbiter. The focus will be on the spectral window containing the
# Oxygen VI line.
# Here we search for a specific observation and query the Solar Orbiter archive for the L2 SPICE data.
res = Fido.search(
    a.Time("2022-04-02 13:00", "2022-04-02 13:25"),
    a.Instrument.spice,
    a.Level(2),
    a.Provider.soar,
    a.soar.Product.spice_n_ras,
)
filename = Fido.fetch(res)[0]

###############################################################################
# Reading the Data
# ----------------
#
# We now open the FITS file directly and pull out the O VI 1032 Å window by its ``EXTNAME``.
# SPICE L2 files store each spectral window as a separate HDU, and "Peak" denotes the narrow window cropped around the line peak.
# This keeps the example simplified, but and shows what's in the file, but for routine SPICE work you'll want the dedicated tools mentioned below.

hdul = fits.open(filename)
hdu = hdul["O VI 1032 - Peak"]

###############################################################################
# If you are using SPICE data more regularly you might wish to look
# into the `read_spice_l2_fits
# <https://docs.sunpy.org/projects/sunraster/en/stable/generated/api/sunraster.instr.spice.read_spice_l2_fits.html>`__
# function and the `sospice
# <https://sospice.readthedocs.io/en/stable/>`__ package for more
# SPICE tools.
#
# The next step is to create an `~ndcube.NDCube` object from the data we have opened.
# To construct this `~ndcube.NDCube` object we use the header from the
# FITS file to make a WCS object, the data and the unit of the data,
# as well as constructing a mask for all points where the data is NaN.

spice = NDCube(
    hdu.data,
    wcs=WCS(hdu),
    unit=hdu.header["BUNIT"],
    mask=np.isnan(hdu.data),
)

###############################################################################
# The first dimension (time / map repeat) is length one so we will drop it.
spice = spice.squeeze()

###############################################################################
# We also crop to a spatial subregion so the example runs quickly.
# The dimensions of the remaining cube are (wavelength, slit, raster step),
# so this slice reduces the size in both the spatial axes.
# If you run this example locally you can comment out this line to run the full raster scan.
spice = spice[:, 200:500, 30:-50]

###############################################################################
# This cube we have constructed has three pixel and four world
# dimensions, the first array dimension is wavelength followed by the
# slit dimension and then the rastering dimension. The world axes are
# wavelength, helioprojective latitude and longitude and time, the time coordinate
# increases along the rastering dimension because each raster step is acquired at a slightly later time as it rasters.

print(spice)

###############################################################################
# To aid our analysis we are going to make two rebinned cubes from
# this data, one summed along the wavelength dimension and one of the
# spectra averaged over all spatial pixels.
#
# We use the shortcut value ``-1`` to `~ndcube.NDCube.rebin` to set
# the number of pixels in a bin equal to the number of pixels along
# that axis, meaning that we put that whole axis into one bin.
#
# We then squeeze the resulting cubes to drop all length one axes.

wl_sum = spice.rebin((-1, 1, 1), operation=np.sum).squeeze()
print(wl_sum)

spatial_mean = spice.rebin((1, -1, -1), operation=np.nanmean).squeeze()
print(spatial_mean)

###############################################################################
# We can now use the `~ndcube.NDCube` built-in plotting routines to show the
# spatially averaged spectra, which uses `~astropy.visualization.wcsaxes.WCSAxes`.

plt.figure()
ax = spatial_mean.plot(axes_units=[u.nm])
ax.coords[0].set_major_formatter("x.xx")

###############################################################################
# Initial Model
# -------------
#
# We now create an initial model for the O VI line, and then fit the
# average spectra to get a strong initial model guess as input for the per-pixel parallel
# fitting.
# The O VI window the spectrum is well-described by a single emission line on top of a roughly flat background, so our model is a constant continuum (`~astropy.modeling.functional_models.Const1D`) plus one Gaussian (`~astropy.modeling.functional_models.Gaussian1D`) centered near 103.2 nm.

OVI_wave = 103.2 * u.nm

initial_model = (
    m.Const1D(amplitude=5 * spice.unit) +
    m.Gaussian1D(amplitude=40 * spice.unit, mean=OVI_wave, stddev=0.05 * u.nm)
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

fitter = TRFLSQFitter()
average_fit = fitter(
    initial_model,
    spatial_mean.axis_world_coords("em.wl")[0].to(u.nm),
    spatial_mean.data * spatial_mean.unit,
    filter_non_finite=True,
)
print(average_fit)

###############################################################################
# Note that we set the ``filter_non_finite`` keyword to the
# `astropy.modeling.fitting.TRFLSQFitter.__call__` function to tell
# the fitter to ignore any ``NaN`` values rather than error.

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
# * A data array. This can be a numpy array or a dask array, or a
#   `~astropy.nddata.NDData` / `~ndcube.NDCube` object. If it's one of
#   these objects then the data, wcs, mask, data_unit and uncertainty
#   are all extracted from the object and used in place of their
#   respective keyword arguments.
# * A model to fit.
# * A fitter instance.
# * The fitting axis (or axes).
#
# What is returned from `~astropy.modeling.fitting.parallel_fit_dask`
# is a model with array parameters with the shape of the non-fitting
# axes of the data, so in this case the array size of the non-wavelength axes.

###############################################################################
# We can therefore fit all our SPICE cube as follows, note we are
# still passing the ``filter_non_finite`` argument as before.

spice_model_fit = parallel_fit_dask(
    data=spice,
    model=average_fit,
    fitter=TRFLSQFitter(),
    fitting_axes=0,
    fitter_kwargs={"filter_non_finite": True},  # Filter out non-finite values,
)

###############################################################################
# We now create a nice stacked plot of the original data and the shift
# in the peak locations of the Gaussian fit. We shall talk more about
# this later.

g1_peak_shift = spice_model_fit.mean_1.quantity.to(
    u.km / u.s,
    equivalencies=u.doppler_optical(OVI_wave),
)

# SPICE pixels are non-square (along-slit ``CDELT2``< raster step
# ``CDELT1``), so we set the aspect accordingly.  We can calculate the
# aspect ratio from the physical size of the pixels to scale the plot
# to the world coordinates.
aspect = hdu.header["CDELT2"] / hdu.header["CDELT1"]

# Create a figure with one column and two rows using the WCS for the
# ``wl_sum`` cube.  Use the constrained layout for better use of space
# in the figure.
fig1, ax1 = plt.subplots(
    subplot_kw=dict(projection=wl_sum),
    figsize=(4.5, 4),
    layout="constrained",
)
fig1.suptitle(f"SPICE - {hdu.header['EXTNAME']} - {hdu.header['DATE-AVG']}")

wl_sum.plot(axes=ax1, norm=LogNorm(), aspect=aspect)
fig1.colorbar(
    ax1.get_images()[0],
    ax=ax1,
    extend="both",
    label=f"{wl_sum.unit:latex}",
    shrink=0.9,
)
ax1.set_title("Data (summed over wavelength)", pad=40)

fig2, ax2 = plt.subplots(
    subplot_kw=dict(projection=wl_sum),
    figsize=(4.5, 4),
    layout="constrained",
)
fig2.suptitle(f"SPICE - {hdu.header['EXTNAME']} - {hdu.header['DATE-AVG']}")

g1_max = np.nanpercentile(np.abs(g1_peak_shift.value), 97)
mean_1 = ax2.imshow(
    g1_peak_shift.value,
    cmap="coolwarm",
    norm=CenteredNorm(halfrange=g1_max),
    aspect=aspect,
)
fig2.colorbar(
    mean_1,
    ax=ax2,
    extend="both",
    label=f"Velocity from Doppler shift [{g1_peak_shift.unit:latex}]",
    shrink=0.9,
)
ax2.set_title(f"O VI ({OVI_wave:latex})", pad=40)

for ax in [ax1, ax2]:
    ax.coords[0].set_ticklabel(exclude_overlapping=True)
    ax.coords[0].set_axislabel("Helioprojective Longitude")
    ax.coords[1].set_axislabel("Helioprojective Latitude")
    ax.coords[2].set_ticklabel(exclude_overlapping=True)

####################################################################
# Note that the velocity map above demonstrates the fitting machinery
# but should not be read as science-grade Doppler measurements.  SPICE
# has a known elliptical and rotated PSF :cite:p:`fludra_spice_2021` and
# channel-dependent wavelength calibration offsets that bias raw line
# centroids; see :cite:p:`plowman_spice_psf_2023,plowman_spice_doppler_2026` for correction methods.

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
# we have two models (0 and 1) added together.

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

spice_model_fit.mean_1.quantity.to(
    u.km/u.s,
    equivalencies=u.doppler_optical(OVI_wave),
)

###############################################################################
# One other thing we may want to do is to evaluate the fitted model for all pixels,
# for example, to plot them or otherwise inspect a single fit.
# This can be done by passing in a wavelength array which is
# :external+numpy:doc:`broadcastable <user/basics.broadcasting>` to
# the shape of the non-fitting axes. We can do this by once again
# using the `~ndcube.NDCube.axis_world_coords` method of
# `~ndcube.NDCube`.

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
# Finally we can make a plot of every per-pixel fit  on top of the spatial-average spectrum and its model.
# The spread of red curves around the dashed line is a visual summary of how the line shifts, broadens,
# and changes intensity across the FOV.
# Notice that it looks like a few pixels have failed to converge on a physical fit.
# These pixels could be excluded, or other things such as
# :external+astropy:doc:`constraints
# <modeling/example-fitting-constraints>`, could be included in the
# fitting to control the fit.

fig = plt.figure(figsize=(11, 5))

ax = spatial_mean.plot(axes_units=[u.nm], label="Average spectra", zorder=99)
ax.coords["wavelength"].set_major_formatter("x.xx")

ax.plot(
    average_fit(spatial_mean.axis_world_coords("em.wl")[0]),
    linestyle="--",
    label="Spatial Average Model",
    zorder=99,
)

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
