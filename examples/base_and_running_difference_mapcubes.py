"""
==================================================
Base Difference and Running Difference in Mapcubes
==================================================

This example illustrates how to do base and running differencing with a MapCube.
Base differencing uses a fixed reference point when compared to running difference.
"""

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

import sunpy.map
import sunpy.physics.differential_rotation as diffrot
from sunpy.data.sample import AIA_193_CUTOUT01_IMAGE, AIA_193_CUTOUT02_IMAGE, AIA_193_CUTOUT03_IMAGE

###############################################################################
# We create the MapCube using the AIA_193_CUTOUT sample data.
# To create a MapCube, we can call Map directly but add in a keyword to output a MapCube instead.
aiamapcube = sunpy.map.Map(AIA_193_CUTOUT01_IMAGE, AIA_193_CUTOUT02_IMAGE,
                           AIA_193_CUTOUT03_IMAGE, cube=True)

############################################################################
# In case of running difference, we loop through all the maps in the aiamapcube and differentially rotate each map
# in the mapcube with respect to the previous map, while in case of base difference we only differentially 
# rotate each map in the mapcube to the time of the base map.
# We then store all such difference maps in a list.
base_diffmap = []
running_diffmap = []
for i, map_i in enumerate(aiamapcube[1:]):
    aiamap_rot = diffrot.diffrot_map(map_i, time=aiamapcube[0].date)
    aiamapcube_rot = diffrot.diffrot_map(aiamapcube[i+1], time=aiamapcube[i].date)
    diffdata = map_i.data - aiamap_rot.data
    smap_base = sunpy.map.Map(diffdata, map_i.meta)
    diffdata = aiamapcube_rot.data - map_i.data
    smap_run = sunpy.map.Map(diffdata, map_i.meta)
    smap_base.plot_settings['cmap'] = cm.get_cmap('Greys_r')
    smap_base.plot_settings['norm'] = colors.LogNorm(100, smap.max())
    smap_run.plot_settings['cmap'] = cm.get_cmap('Greys_r')
    smap_run.plot_settings['norm'] = colors.LogNorm(100, smap.max())
    base_diffmap.append(smap_base)
    running_diffmap.append(smap_run)

############################################################################
# This plots the original mapcube
aiamapcube.peek()

##############################################################################
# This plots the final mapcube obtained after implementing the base difference
result_mapcube = sunpy.map.MapCube(base_diffmap)
result_mapcube.peek()

############################################################################
# This plots the final mapcube after implementing the running difference
result_mapcube = sunpy.map.MapCube(running_diffmap)
result_mapcube.peek()
plt.show()
