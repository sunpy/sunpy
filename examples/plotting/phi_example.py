"""
===============================
Plotting Solar Orbiter PHI Data
===============================

This example demonstrates how to plot Solar Orbiter PHI data.
"""

import matplotlib.pyplot as plt

from astropy.time import Time

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# Searching for PHI-HRT Data
# --------------------------
#
# We first search for all **Solar Orbiter PHI-HRT** (High Resolution Telescope) data products # (excluding 'phi-hrt-stokes' which is not compatible with sunpy.map)

t_start_hrt = Time('2024-03-23T20:00', format='isot', scale='utc')
t_end_hrt = Time('2024-03-23T23:59', format='isot', scale='utc')

search_results_phi_hrt_all = Fido.search(a.Instrument('PHI'), a.Time(t_start_hrt.value, t_end_hrt.value), (a.soar.Product('phi-hrt-blos') | a.soar.Product('phi-hrt-bmag') | a.soar.Product('phi-hrt-binc') | a.soar.Product('phi-hrt-bazi')| a.soar.Product('phi-hrt-vlos') | a.soar.Product('phi-hrt-icnt')))

print(search_results_phi_hrt_all)

###############################################################################
# Fetching the First Available PHI-HRT Files
# ------------------------------------------
files_phi_hrt_all = Fido.fetch(search_results_phi_hrt_all[:, 0])


###############################################################################
# Loading and Plotting PHI-HRT Data
# ---------------------------------

phi_hrt_blos_map = sunpy.map.Map(files_phi_hrt_all[0])
phi_hrt_bmag_map = sunpy.map.Map(files_phi_hrt_all[1])
phi_hrt_binc_map = sunpy.map.Map(files_phi_hrt_all[2])
phi_hrt_bazi_map = sunpy.map.Map(files_phi_hrt_all[3])
phi_hrt_vlos_map = sunpy.map.Map(files_phi_hrt_all[4])
phi_hrt_icnt_map = sunpy.map.Map(files_phi_hrt_all[5])

########################################

phi_hrt_blos_map.plot()
plt.colorbar(label=phi_hrt_blos_map.unit.to_string())

########################################

phi_hrt_bmag_map.plot()
plt.colorbar(label=phi_hrt_bmag_map.unit.to_string())
plt.clim(0,2500)

########################################

phi_hrt_binc_map.plot()
plt.colorbar(label=phi_hrt_binc_map.unit.to_string())

########################################

phi_hrt_bazi_map.plot()
plt.colorbar(label=phi_hrt_bazi_map.unit.to_string())

########################################

phi_hrt_vlos_map.plot()
plt.colorbar(label=phi_hrt_vlos_map.unit.to_string())
plt.clim(-2,2)

########################################

phi_hrt_icnt_map.plot()
plt.colorbar(label=phi_hrt_icnt_map.unit.to_string())
plt.clim(0,1.2)

###############################################################################
# Searching for PHI-FDT Data
# --------------------------------
#
# We first search for all **Solar Orbiter PHI-FDT** (Full Disc Telescope) data products
# For this time range, at the time of this example creation 'phi-fdt-vlos' is not available yet via SOAR
# (excluding 'phi-fdt-stokes' which is not compatible with sunpy.map)

t_start_fdt = Time('2025-02-25T20:00', format='isot', scale='utc')
t_end_fdt = Time('2025-02-25T23:59', format='isot', scale='utc')

search_results_phi_fdt_all = Fido.search(a.Instrument('PHI'), a.Time(t_start_fdt.value, t_end_fdt.value), (a.soar.Product('phi-fdt-blos') | a.soar.Product('phi-fdt-bmag') | a.soar.Product('phi-fdt-binc') | a.soar.Product('phi-fdt-bazi')| a.soar.Product('phi-fdt-icnt')))

print(search_results_phi_fdt_all)

###############################################################################
# Fetching the First Available PHI-fdt File
# -----------------------------------------
files_phi_fdt_all = Fido.fetch(search_results_phi_fdt_all)

###############################################################################
# Loading and Plotting PHI-FDT Data
# ---------------------------------
# We load each map and rotate to point Solar North up
phi_fdt_blos_map = sunpy.map.Map(files_phi_fdt_all[0]).rotate(recenter=True)
phi_fdt_bmag_map = sunpy.map.Map(files_phi_fdt_all[1]).rotate(recenter=True)
phi_fdt_binc_map = sunpy.map.Map(files_phi_fdt_all[2]).rotate(recenter=True)
phi_fdt_bazi_map = sunpy.map.Map(files_phi_fdt_all[3]).rotate(recenter=True)
phi_fdt_icnt_map = sunpy.map.Map(files_phi_fdt_all[4]).rotate(recenter=True)

###############################################################################
# Make a mask to mask out off-disc pixels and make new sunpy maps
# -----------------------------------------
hpc_coords = sunpy.map.all_coordinates_from_map(phi_fdt_blos_map)
mask = ~sunpy.map.coordinate_is_on_solar_disk(hpc_coords)

phi_fdt_blos_map = sunpy.map.Map(phi_fdt_blos_map.data, phi_fdt_blos_map.meta, mask=mask)
phi_fdt_bmag_map = sunpy.map.Map(phi_fdt_bmag_map.data, phi_fdt_bmag_map.meta, mask=mask)
phi_fdt_binc_map = sunpy.map.Map(phi_fdt_binc_map.data, phi_fdt_binc_map.meta, mask=mask)
phi_fdt_bazi_map = sunpy.map.Map(phi_fdt_bazi_map.data, phi_fdt_bazi_map.meta, mask=mask)
phi_fdt_icnt_map = sunpy.map.Map(phi_fdt_icnt_map.data, phi_fdt_icnt_map.meta, mask=mask)

########################################

phi_fdt_blos_map.plot()
plt.colorbar(label=phi_fdt_blos_map.unit.to_string())

########################################

phi_fdt_bmag_map.plot()
plt.colorbar(label=phi_fdt_bmag_map.unit.to_string())
plt.clim(0,1500)

########################################

phi_fdt_binc_map.plot()
plt.colorbar(label=phi_fdt_binc_map.unit.to_string())

########################################

phi_fdt_bazi_map.plot()
plt.colorbar(label=phi_fdt_bazi_map.unit.to_string())

########################################

phi_fdt_icnt_map.plot()
plt.colorbar(label=phi_fdt_icnt_map.unit.to_string())
plt.clim(0,1.2)
