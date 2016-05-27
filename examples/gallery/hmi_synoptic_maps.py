"""
HMI Daily Synoptic Map
----------------------

In this example we will load a FITS file not directly supported by SunPy, and
build a custom coordinate frame to display the coordinate system correctly.

This notebook plots the HMI Daily Synoptic Maps, the file used in this example
can be downloaded from
[here](http://jsoc.stanford.edu/data/hmi/synoptic/hmi.Mldailysynframe_720s_nrt.fits).
"""

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
import astropy.coordinates as coord
from astropy.io import fits
from astropy.coordinates import frame_transform_graph
from astropy.utils.data import download_file

import sunpy.map
import sunpy.coordinates as spc

###############################################################################
# Use astropy to download the file to a temp location.

filename = download_file('http://jsoc.stanford.edu/data/hmi/synoptic/hmi.Mldailysynframe_720s_nrt.fits')


###############################################################################
# We read the file in using astropy to give the opertunity to fix it.

hdus = fits.open(filename)
hdus.verify('fix')

# extract the data
data = hdus[0].data


###############################################################################
# Convert the header to a dict to make modifcation easier.

header = dict(hdus[0].header)


###############################################################################
# There are a couple of oddities with this file, firstly the value of 'CUNIT2':

print(header['CUNIT2'])


###############################################################################
# That is not a unit! What this is telling us is that the latitude coordinate
# is actually the sine of latitude, this is going to cause us some issues
# later, but first let us fix that so SunPy can read it:

header['CUNIT2'] = 'deg'


###############################################################################
# Also the value of 'HGLN_OBS':

print(header['HGLN_OBS'])


###############################################################################
# Is not a value SunPy will parse, so we delete it and SunPy will default it to
# zero:

del header['HGLN_OBS']


###############################################################################
# Now we create a SunPy Map from the data and header:

m = sunpy.map.Map((data, header))
# Set the colorbar properties.
m.plot_settings['cmap'] = 'hmimag'
m.plot_settings['norm'] = plt.Normalize(-1500, 1500)


###############################################################################
# Now we have to consider how we are going to represent the values of solar
# latitude on the graph when the image is in sine latitude. First, let's
# inspect the coordinate system of the image:

print(m.coordinate_system)


###############################################################################
# This is Carrington Heliographic with the Cylindircal Equal Area projection,
# which SunPy will recognise:

print(m.coordinate_frame)


###############################################################################
# This poses a problem, as the coordinate frame is not actually Heliographic
# Carrington as the latitude is actually sine latitude. However, SunPy will
# happily represent this data ok, but if we want to display the real latitude
# we will need to define a conversion from what SunPy interprets to be
# Heliographic Carrington to this Sine HGC frame. To do this we define a
# `ASineHGC` frame, this is the reverse transformation, i.e. SunPy thinks it is
# from HGC to ASineHGC but is actually from SineHGC to HGC. We do this by
# creating a subclass of the `HelioGraphicCarrington` frame, which has to
# change nothing (the units of sine latitude are faked to be in degrees or else
# it all breaks) and then define two transformation functions.


class ASineHGC(spc.HeliographicCarrington):
    """
    This frame is a reverse transformation from the HGC frame with CUNIT2
    sin(deg) to the real HGC frame. It is defined backwards because the WCS
    object for the HGC sin(deg) frame reads in as a proper HGC frame due to the
    values of CTYPE. This frame therefore converts from what WCS thinks is HGC
    to ASineHGC but is in reality converting from ASineHGC to HGC.
    """


@frame_transform_graph.transform(coord.FunctionTransform,
                                 spc.HeliographicCarrington, ASineHGC)
def hgc_to_sinehgc(hgc_coord, sinehgc_frame):
    lat = hgc_coord.lat
    lon = hgc_coord.lon

    lat_out = u.Quantity(np.arcsin(lat.value), u.rad)

    return sinehgc_frame.realize_frame(
        coord.UnitSphericalRepresentation(lat=lat_out, lon=lon))


@frame_transform_graph.transform(coord.FunctionTransform, ASineHGC,
                                 spc.HeliographicCarrington)
def sinehgc_to_hgc(sinehgc_coord, hgc_frame):
    lat = sinehgc_coord.lat
    lon = sinehgc_coord.lon

    lat_out = u.Quantity(np.sin(lat.value), u.deg)

    return hgc_frame.realize_frame(
        coord.UnitSphericalRepresentation(lat=lat_out, lon=lon))


###############################################################################
# Once we have this defined we can create a plot using SunPy and WCSAxes:


###############################################################################
# Create a figure with the Map's projection:
fig = plt.figure(figsize=(12, 5))
axes = plt.subplot(projection=m)

# Plot the image
im = m.plot()

# Set up the Sine Latitude Grid
x = axes.coords[0]
y = axes.coords[1]

x.set_coord_type('longitude', coord_wrap=360.)

x.set_major_formatter('dd')
y.set_major_formatter('d.d')

x.set_axislabel("Carrington Longitude [deg]")
y.set_axislabel("Sine Latitude")


x.set_ticks(color='black', exclude_overlapping=True)
y.set_ticks(color='black', exclude_overlapping=True)

# Hide the grid
axes.coords.grid(alpha=0)

# Create a colorbar
cb = plt.colorbar(im, fraction=0.019, pad=0.1)
cb.set_label("LOS Magnetic Field [Gauss]")

# Now create the overlay in the actual HGC coordinate frame

overlay = axes.get_coords_overlay('asinehgc')

lon = overlay[0]
lat = overlay[1]

lon.set_ticklabel_visible(False)

lat.set_major_formatter('dd')

lat.set_axislabel('Solar Latitude [deg]')

lat.set_ticks_position('tr')
lat.set_ticks(spacing=10*u.deg, exclude_overlapping=True)

# Another horrible hack to make the ticks draw on the RHS
axes.set_xlim((0, 3585))

plt.title("HMI Daily Synoptic Frame for Carrington Rotation"
          " {}-{}".format(header['CAR_ROT'], header['CAR_ROT']+1))

plt.show()
