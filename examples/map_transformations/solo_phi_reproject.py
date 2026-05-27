"""
================================================================
Downloading and Reprojecting a Solar Orbiter PHI/HRT Magnetogram
================================================================

This example demonstrates querying for Solar Orbiter PHI data and reprojecting it.
"""

# sphinx_gallery_tags = ["Acquiring Data", "Solar Orbiter",  "SOAR"]

import matplotlib.pyplot as plt
from matplotlib.colors import CenteredNorm

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.coordinates
import sunpy.map
import sunpy.visualization.colormaps
from sunpy.map.header_helper import make_fitswcs_header
from sunpy.net import Fido
from sunpy.net import attrs as a

################################################################################
# Do a search from the SOAR for HRT line of sight magnetic field.
search_results_phi_hrt = Fido.search(
    a.Time("2024-10-14T00:25:00", "2024-10-14T00:35:00"),
    a.soar.Product("phi-hrt-blos"),
)
print(search_results_phi_hrt)

################################################################################
# Download just the first result from the SOAR client.
sr_phi_hrt = search_results_phi_hrt["soar", 0]
blos_file = Fido.fetch(sr_phi_hrt)

################################################################################
# Create a Map and plot with a clipped color bar.
blos_map = sunpy.map.Map(blos_file[0])
blos_map.plot(norm=CenteredNorm(halfrange=1500, clip=True))
plt.colorbar()

################################################################################
# Construct a Heliographic Stonyhurst header to reproject the image to.
# We choose to build a map where the reference coordinate is at grid
# center so that we have a rectilinear grid.

hgs_center = SkyCoord(0, 0, unit=u.arcsec, frame="heliographic_stonyhurst", obstime=blos_map.reference_date)
cea_scale = (0.03, 0.03) * u.deg / u.pix

lon_width = 50 * u.deg
lat_height = 25 * u.deg

nx = int((lon_width / cea_scale[0]).to_value(u.pix))
ny = int((lat_height / cea_scale[1]).to_value(u.pix))

cea_hdr = make_fitswcs_header(
    (ny, nx),
    hgs_center,
    reference_pixel=(nx / 2 - 200, ny) * u.pix,  # Offset the reference pixel to re-center the image
    scale=cea_scale,
    projection_code="CEA",
    instrument=blos_map.instrument,
    wavelength=blos_map.wavelength,
)

################################################################################
# Reproject the map to the new header
# Note that this treats the LOS magnetic field values as scalar image data.
# This is useful for visualisation, but it does not convert B_LOS into a radial
# or local heliographic magnetic field component.

outmap = blos_map.reproject_to(cea_hdr, algorithm="adaptive", kernel="Hann")

fig = plt.figure()
outmap.plot(norm=CenteredNorm(halfrange=1500, clip=True))

plt.show()
