"""
=====================================
Plotting Solar Orbiter Fields of View
=====================================

In this example we shall demonstrate how to plot the fields of view of multiple Solar Orbiter instruments.
"""
# sphinx_gallery_tags = ["Coordinates", "Solar Orbiter", "Map", "SOAR", "EUI", "SPICE", "AIA"]
# sphinx_gallery_thumbnail_number = 2

import matplotlib.pyplot as plt
import numpy as np
from ndcube import NDCube

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.visualization import drawing

################################################################################
# Download and Read Multiple Observations
# ---------------------------------------
#
# The first step is to search for and download multiple datasets.
# During this time window a coordinated campaign between Solar Orbiter and DKIST (and PSP??!) was underway.
time_range = a.Time("2022-10-24T18:55", "2022-10-24T19:35")

################################################################################
# We search for EUI Full Sun Imager, EUI High Resolution Imager and SPICE data.
solo = (
    a.soar.Product("EUI-HRIEUV174-IMAGE") | a.soar.Product("EUI-FSI174-IMAGE") | a.soar.Product("SPICE-N-RAS")
) & a.Level(2)

################################################################################
# As well as AIA 17.1 nm data.
aia = a.Instrument.aia & a.Wavelength(171 * u.Angstrom) & a.Sample(10 * u.minute)

results = Fido.search(
    time_range,
    solo | aia,
)
print(results[:, 0])

################################################################################
# We will then just download the first results for all our different data sources.
files = Fido.fetch(results[:, 0], site="NSO")
# Sort the files by filename
files.sort()

################################################################################
# Open the first three into sunpy maps
aia, eui_fsi, eui_hri = sunpy.map.Map(files[:-1])

################################################################################
# Load the SPICE file into a NDCube
hdul = fits.open(files[-1])
hdu = hdul["N IV 765 - SH - Comp 8 ... Ne VIII 770 - LH - Comp 8 (Merged)"]

spice = NDCube(
    hdu.data,
    wcs=WCS(hdu),
    unit=hdu.header["BUNIT"],
    mask=np.isnan(hdu.data),
)

# Drop the length one time / raster repeat dimension
spice = spice.squeeze()

spice_wl_sum = spice.rebin((-1, 1, 1), operation=np.sum).squeeze()

################################################################################
# Initial Plots

# fig = plt.figure()
# ax = fig.add_subplot(projection=eui_fsi)
# eui_fsi.plot(axes=ax)

################################################################################
# We can see that the FSI image is a little zoomed out, so let's crop
# it down to closer to the limb.

eui_fsi_zoom = eui_fsi.submap(
    bottom_left=SkyCoord(-3000, -3000, unit="arcsec", frame=eui_fsi.coordinate_frame),
    top_right=SkyCoord(3000, 3000, unit="arcsec", frame=eui_fsi.coordinate_frame),
)

################################################################################
# We shall now plot the extent of EUI HRI and SPICE, as well as the limb as seen from AIA on this image.

fig = plt.figure()
ax = fig.add_subplot(projection=eui_fsi_zoom)
eui_fsi_zoom.plot(axes=ax)
eui_hri.draw_extent(label="EUI HRI")
drawing.extent(axes=ax, wcs=spice_wl_sum.wcs, color="blue", label="SPICE")
drawing.limb(ax, aia.observer_coordinate, rsun=aia.rsun_meters, color='C3', lw=2, label='AIA limb')
ax.legend()

################################################################################
# We shall now do the same but based on the AIA image.

fig = plt.figure()
ax = fig.add_subplot(projection=aia)
aia.plot(axes=ax)
eui_hri.draw_extent(label="EUI HRI")
drawing.extent(axes=ax, wcs=spice_wl_sum.wcs, color="blue", label="SPICE")
drawing.limb(ax, eui_fsi.observer_coordinate, rsun=eui_fsi.rsun_meters, color='C3', lw=2, label='EUI limb')
ax.legend()

plt.show()
