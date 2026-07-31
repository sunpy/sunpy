"""
========================================
Combining Datasets and Plotting PHI Data
========================================

This example shows how to combine various PHI datasets and plot them together.
The example uses data from the PHI instrument on board Solar Orbiter.

.. warning::

    This example requires the `sunkit-image <https://docs.sunpy.org/projects/sunkit-image/en/latest/>`__ package.
"""
# sphinx_gallery_tags = ["Solar Orbiter", "PHI", "Visualization"]

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.map
from sunpy.coordinates.utils import coordinate_is_on_solar_disk
from sunpy.data.sample import (
    PHI_FDT_BINC_IMAGE,
    PHI_FDT_BAZI_IMAGE,
    PHI_FDT_BMAG_IMAGE,
    PHI_FDT_ICNT_IMAGE,
    PHI_FDT_VLOS_IMAGE,
)
from sunpy.map.maputils import all_coordinates_from_map

###############################################################################
# We start by loading each image. The individual maps can be downloaded using
# the sample data downloader.

phi_fdt_bmag_map = sunpy.map.Map(PHI_FDT_BMAG_IMAGE)
phi_fdt_binc_map = sunpy.map.Map(PHI_FDT_BINC_IMAGE)
phi_fdt_bazi_map = sunpy.map.Map(PHI_FDT_BAZI_IMAGE)
phi_fdt_vlos_map = sunpy.map.Map(PHI_FDT_VLOS_IMAGE)
phi_fdt_icnt_map = sunpy.map.Map(PHI_FDT_ICNT_IMAGE)

###############################################################################
# The data contains NaN values off the solar disk. We now create a mask using
# the ``coordinate_is_on_solar_disk`` function. This is a useful function for
# creating masks, as we have done in other examples such as
# :ref:`sphx_glr_generated_gallery_plotting_masked_composite_plot.py`.

hpc_coords = all_coordinates_from_map(phi_fdt_bmag_map)
mask = ~coordinate_is_on_solar_disk(hpc_coords)

###############################################################################
# We now apply the mask to each of the maps.

phi_fdt_bmag_map = sunpy.map.Map(phi_fdt_bmag_map.data, phi_fdt_bmag_map.meta, mask=mask)
phi_fdt_binc_map = sunpy.map.Map(phi_fdt_binc_map.data, phi_fdt_binc_map.meta, mask=mask)
phi_fdt_bazi_map = sunpy.map.Map(phi_fdt_bazi_map.data, phi_fdt_bazi_map.meta, mask=mask)
phi_fdt_vlos_map = sunpy.map.Map(phi_fdt_vlos_map.data, phi_fdt_vlos_map.meta, mask=mask)
phi_fdt_icnt_map = sunpy.map.Map(phi_fdt_icnt_map.data, phi_fdt_icnt_map.meta, mask=mask)

###############################################################################
# Next, we plot each of the maps using a colormap appropriate for the
# respective quantity. We also add a colorbar to each plot.

fig = plt.figure(figsize=(10, 15))

########################################

norm = ImageNormalize(vmin=0, vmax=1500, stretch=SqrtStretch(), clip=False)
ax1 = fig.add_subplot(5, 1, 1, projection=phi_fdt_bmag_map)
transparent_cmap = plt.get_cmap("cubehelix")
transparent_cmap.set_bad(alpha=0.0)
phi_fdt_bmag_map.plot(axes=ax1, cmap=transparent_cmap, norm=norm, annotate=False)
ax1.set_title(f"{phi_fdt_bmag_map.instrument} {phi_fdt_bmag_map.measurement}")
ax1.set_xlabel("Solar-X [arcsec]")
ax1.set_ylabel("Solar-Y [arcsec]")
plt.colorbar(label=phi_fdt_bmag_map.unit.to_string())
plt.clim(0, 1500)

########################################

ax2 = fig.add_subplot(5, 1, 2, projection=phi_fdt_binc_map)
phi_fdt_binc_map.plot(axes=ax2)
ax2.set_title(f"{phi_fdt_binc_map.instrument} {phi_fdt_binc_map.measurement}")
ax2.set_xlabel("Solar-X [arcsec]")
ax2.set_ylabel("Solar-Y [arcsec]")
plt.colorbar(label=phi_fdt_binc_map.unit.to_string())

########################################

ax3 = fig.add_subplot(5, 1, 3, projection=phi_fdt_bazi_map)
phi_fdt_bazi_map.plot(axes=ax3)
ax3.set_title(f"{phi_fdt_bazi_map.instrument} {phi_fdt_bazi_map.measurement}")
ax3.set_xlabel("Solar-X [arcsec]")
ax3.set_ylabel("Solar-Y [arcsec]")
plt.colorbar(label=phi_fdt_bazi_map.unit.to_string())

########################################

ax4 = fig.add_subplot(5, 1, 4, projection=phi_fdt_vlos_map)
phi_fdt_vlos_map.plot(axes=ax4)
ax4.set_title(f"{phi_fdt_vlos_map.instrument} {phi_fdt_vlos_map.measurement}")
ax4.set_xlabel("Solar-X [arcsec]")
ax4.set_ylabel("Solar-Y [arcsec]")
plt.colorbar(label=phi_fdt_vlos_map.unit.to_string())
plt.clim(-2, 2)

########################################

ax5 = fig.add_subplot(5, 1, 5, projection=phi_fdt_icnt_map)
phi_fdt_icnt_map.plot(axes=ax5)
ax5.set_title(f"{phi_fdt_icnt_map.instrument} {phi_fdt_icnt_map.measurement}")
ax5.set_xlabel("Solar-X [arcsec]")
ax5.set_ylabel("Solar-Y [arcsec]")
plt.colorbar(label=phi_fdt_icnt_map.unit.to_string())
plt.clim(0, 1.2)

plt.show()
