"""
===============================
Plotting Solar Orbiter PHI Data
===============================

This example demonstrates how to plot Solar Orbiter PHI data.
"""

import matplotlib.pyplot as plt

from astropy.time import Time

import sunpy.map
import sunpy.visualization.colormaps
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# Searching for PHI-HRT Data
# --------------------------------
#
# We first search for all **Solar Orbiter PHI-HRT** (High Resolution Telescope) data products

t_start_hrt = Time('2024-03-23T20:00', format='isot', scale='utc')
t_end_hrt = Time('2024-03-23T23:59', format='isot', scale='utc')

search_results_phi_hrt_all = Fido.search(a.Instrument('PHI'), a.Time(t_start_hrt.value, t_end_hrt.value), (a.soar.Product('phi-hrt-blos') | a.soar.Product('phi-hrt-bmag') | a.soar.Product('phi-hrt-binc') | a.soar.Product('phi-hrt-bazi')| a.soar.Product('phi-hrt-vlos') | a.soar.Product('phi-hrt-icnt')))

print(search_results_phi_hrt_all)

###############################################################################
# Fetching the First Available PHI-HRT Files
# -----------------------------------------
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
