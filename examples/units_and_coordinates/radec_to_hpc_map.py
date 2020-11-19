"""
==============================================================================
Create a Helioprojective Map from observations in the RA-DEC coordinate system
==============================================================================

How to create a `~sunpy.map.Map` in Helioprojective Coordinate Frame from radio observations
in GCRS (RA-DEC).

In this example a LOFAR FITS file (created with LOFAR's `Default Pre-Processing Pipeline (DPPP) and
WSClean Imager <https://support.astron.nl/LOFARImagingCookbook/dppp.html>`__) is read in,
the WCS header information is then used to make a new header with the information in Helioprojective,
and a `~sunpy.map.Map` is made.

The LOFAR example file has a WCS in celestial coordinates i.e. Right Ascension and
Declination (RA-DEC). For this example, we are assuming that the definition of LOFAR's
coordinate system for this observation is exactly the same as Astropy's ~astropy.coordinates.GCRS.
For many solar studies we may want to plot this data in some Sun-centered coordinate frame,
such as `~sunpy.coordinates.frames.Helioprojective`. In this example we read the data and
header information from the LOFAR FITS file and then create a new header with updated WCS
information to create a `~sunpy.map.Map` with a HPC coordinate frame. We will make use of the
`astropy.coordinates` and `sunpy.coordinates` submodules together with `~sunpy.map.make_fitswcs_header`
to create a new header and generate a `~sunpy.map.Map`.
"""
import matplotlib.pyplot as plt
import numpy as np

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time

import sunpy.data.sample
import sunpy.map
from sunpy.coordinates import frames, sun

##############################################################################
# We will first begin be reading in the header and data from the FITS file.
hdu = fits.open(sunpy.data.sample.LOFAR_IMAGE)
header = hdu[0].header

#####################################################################################
# The data in this file is in a datacube structure to hold difference frequencies and
# polarizations. We are only interested in the image at one frequency (the only data
# in the file) so we index the data array to be the 2D data of interest.
data = hdu[0].data[0, 0, :, :]

################################################################################
# We can inspect the header, for example, we can print the coordinate system and
# type projection, which here is RA-DEC.
print(header['ctype1'], header['ctype2'])

###############################################################################
# Lets pull out the observation time and wavelength from the header, we will use
# these to create our new header.
obstime = Time(header['date-obs'])
frequency = header['crval3']*u.Hz


###############################################################################
# To create a new `~sunpy.map.Map` header we need convert the reference coordinate
# in RA-DEC (that is in the header) to Helioprojective. To do this we will first create
# an `astropy.coordinates.SkyCoord` of the reference coordinate from the header information.
# We will need the location of the observer (i.e. where the observation was taken).
# We first establish the location on Earth from which the observation takes place, in this
# case LOFAR observations are taken from Exloo in the Netherlands, which we define in lat and lon.
# We can convert this to a SkyCoord in GCRSat the observation time.

lofar_loc = EarthLocation(lat=52.905329712*u.deg, lon=6.867996528*u.deg)
lofar_gcrs = SkyCoord(lofar_loc.get_gcrs(obstime))

##########################################################################
# We can then define the reference coordinate in terms of RA-DEC from the header information.
# Here we are using the `obsgeoloc` keyword argument to take into account that the observer is not
# at the center of the Earth (i.e. the GCRS origin). The distance here is the Sun-observer distance.
reference_coord = SkyCoord(header['crval1']*u.Unit(header['cunit1']),
                           header['crval2']*u.Unit(header['cunit2']),
                           frame='gcrs',
                           obstime=obstime,
                           obsgeoloc=lofar_gcrs.cartesian,
                           obsgeovel=lofar_gcrs.velocity.to_cartesian(),
                           distance=lofar_gcrs.hcrs.distance)

##########################################################################
# Now we can convert the `reference_coord` to the HPC coordinate frame
reference_coord_arcsec = reference_coord.transform_to(frames.Helioprojective(observer=lofar_gcrs))

##########################################################################
# Now we need to get the other parameters from the header that will be used
# to create the new header - here we can get the cdelt1 and cdelt2 which are
# the spatial scales of the data axes.
cdelt1 = (np.abs(header['cdelt1'])*u.deg).to(u.arcsec)
cdelt2 = (np.abs(header['cdelt2'])*u.deg).to(u.arcsec)

##################################################################################
# Finally, we need to specify the orientation of the HPC coordinate grid because
# GCRS north is not in the same direction as HPC north.  For convenience, we use
# :func:`~sunpy.coordinates.sun.P` to calculate this relative rotation angle,
# although due to subtleties in definitions, the returned value is inaccurate by
# 2 arcmin, equivalent to a worst-case shift of 0.6 arcsec for HPC coordinates
# on the disk. The image will need to be rotated by this angle so that solar north
# is pointing up.
P1 = sun.P(obstime)

##########################################################################
# Now we can use this information to create a new header using the helper
# function `~sunpy.map.make_fitswcs_header()`. This will create a MetaDict
# which we contain all the necessay WCS information to create a `~sunpy.map.Map`.
# We provide a reference coordinate (in HPC), the spatial
# scale of the observation (i.e. `cdelt1` and `cdelt2`), and the rotation angle (P1).
# Note that here, 1 is subtracted from the crpix1 and crpix2 values, this is because
# the `reference_pixel` keyword in ~sunpy.map.make_fitswcs_header` is zero indexed rather
# than the fits convention of 1 indexed.

new_header = sunpy.map.make_fitswcs_header(data, reference_coord_arcsec,
                                           reference_pixel=u.Quantity([header['crpix1']-1,
                                                                       header['crpix2']-1]*u.pixel),
                                           scale=u.Quantity([cdelt1, cdelt2]*u.arcsec/u.pix),
                                           rotation_angle=-P1,
                                           wavelength=frequency.to(u.MHz).round(2),
                                           observatory='LOFAR')

##########################################################################
# Let's inspect the new header
print(new_header)

##########################################################################
# Lets create a `~sunpy.map.Map`
lofar_map = sunpy.map.Map(data, new_header)

##########################################################################
# We can now plot this map and inspect it
fig = plt.figure()
ax = fig.add_subplot(projection=lofar_map)
lofar_map.plot(cmap='viridis')
lofar_map.draw_limb()
plt.show()

##########################################################################
# We can now rotate the image so that solar north is pointing up and create
# a submap in the field of view of interest.
lofar_map_rotate = lofar_map.rotate()
bl = SkyCoord(-1500*u.arcsec, -1500*u.arcsec, frame=lofar_map_rotate.coordinate_frame)
tr = SkyCoord(1500*u.arcsec, 1500*u.arcsec, frame=lofar_map_rotate.coordinate_frame)
lofar_submap = lofar_map_rotate.submap(bl, top_right=tr)

##########################################################################
# Now lets plot this map, and overplot some contours
fig = plt.figure()
ax = fig.add_subplot(projection=lofar_submap)
lofar_submap.plot(cmap='viridis')
lofar_submap.draw_limb()
lofar_submap.draw_grid()
lofar_submap.draw_contours(np.arange(30, 100, 5)*u.percent)
plt.show()
