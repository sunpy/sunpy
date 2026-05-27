"""
=====================================
Plotting Solar Orbiter Fields of View
=====================================

In this example we shall demonstrate how to plot the fields of view of multiple Solar Orbiter instruments.
"""
# sphinx_gallery_tags = ["Coordinates", "Solar Orbiter", "Map", "SOAR", "EUI", "SPICE", "AIA", "DKIST", "VBI", "VISP"]
# sphinx_gallery_thumbnail_number = -1

import dkist.net  # NOQA: F401
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

dkist_visp = a.Instrument.visp

dkist_vbi = a.Instrument.vbi & a.Wavelength(430 * u.nm, 431 * u.nm)

results = Fido.search(
    time_range,
    solo | aia | dkist_vbi | dkist_visp,
)
print(results[:, 0])

################################################################################
# We will then just download the first results for all our different data sources.
files = Fido.fetch(results[:, 0], site="NSO")
# Sort the files by filename, but lowercase
files.sort(key=str.lower)
# Put the file paths into a dict
files = {name: path for name, path in zip(["AIA", "EUI-FSI", "EUI-HRI", "SPICE", "VBI", "VISP"], files)}

################################################################################
# Open the images as sunpy maps
aia, eui_fsi, eui_hri = sunpy.map.Map(files["AIA"], files["EUI-FSI"], files["EUI-HRI"])

################################################################################
# Load the SPICE file into a NDCube
hdul = fits.open(files["SPICE"])
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
# Open the VBI and VISP ASDF files with the dkist package

vbi = dkist.load_dataset(files["VBI"])
visp = dkist.load_dataset(files["VISP"])

visp_I_wl_sum = visp[0].rebin((1, -1, 1), operation=np.sum).squeeze()

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

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection=eui_fsi_zoom)
eui_fsi_zoom.plot(axes=ax)
eui_hri.draw_extent(label="EUI HRI", color="C1")
drawing.extent(axes=ax, wcs=spice_wl_sum.wcs, color="C2", label="SPICE")

# Add the AIA limb
visible, hidden = drawing.limb(ax, aia.observer_coordinate, rsun=aia.rsun_meters, color="C3", lw=2)
visible.set_label("AIA limb")
hidden.set_label("AIA limb (hidden)")

# Add the VISP FOV
visible, hidden = drawing.extent(ax, wcs=visp_I_wl_sum.wcs, color="C4")
visible.set_label("VISP")

# The VBI data is a mosaic of 9 images, plot each one
for ds in vbi.flat:
    visible, hidden = drawing.extent(ax, ds.wcs, color="C5")
# Only add the last one to the legend
visible.set_label("VBI")

ax.legend()

################################################################################
# We shall now do the same but based on the AIA image.

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection=aia)
aia.plot(axes=ax)
eui_hri.draw_extent(label="EUI HRI")
drawing.extent(axes=ax, wcs=spice_wl_sum.wcs, color="C2", label="SPICE")

# Add the EUI FSI limb
visible, hidden = drawing.limb(ax, eui_fsi.observer_coordinate, rsun=eui_fsi.rsun_meters, color="C3", lw=2)
visible.set_label("EUI limb")
hidden.set_label("EUI limb (hidden)")

# Add the VISP FOV
visible, hidden = drawing.extent(ax, wcs=visp_I_wl_sum.wcs, color="C4")
visible.set_label("VISP")

# The VBI data is a mosaic of 9 images, plot each one
for ds in vbi.flat:
    visible, hidden = drawing.extent(ax, ds.wcs, color="C5")
# Only add the last one to the legend
visible.set_label("VBI")

ax.legend()

################################################################################
# And finally let's zoom in on the shared FOV

eui_fsi_crop = eui_fsi.submap(
    bottom_left=SkyCoord(-100, -500, unit="arcsec", frame=eui_fsi.coordinate_frame),
    top_right=SkyCoord(2000, 1500, unit="arcsec", frame=eui_fsi.coordinate_frame),
)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection=eui_fsi_crop)
eui_fsi_crop.plot(axes=ax)
eui_hri.draw_extent(label="EUI HRI", color="C1")
drawing.extent(axes=ax, wcs=spice_wl_sum.wcs, color="C2", label="SPICE")

# Add the AIA limb
visible, hidden = drawing.limb(ax, aia.observer_coordinate, rsun=aia.rsun_meters, color="C3", lw=2)
visible.set_label("AIA limb")
hidden.set_label("AIA limb (hidden)")

# Add the VISP FOV
visible, hidden = drawing.extent(ax, wcs=visp_I_wl_sum.wcs, color="C4")
visible.set_label("VISP")

# The VBI data is a mosaic of 9 images, plot each one
for ds in vbi.flat:
    visible, hidden = drawing.extent(ax, ds.wcs, color="C5")
# Only add the last one to the legend
visible.set_label("VBI")

ax.legend()

plt.show()
