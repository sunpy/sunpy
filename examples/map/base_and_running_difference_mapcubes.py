"""
===========================================
Base and Running Difference in MapSequences
===========================================

This example illustrates how to do base and running differencing with a MapSequence.
Base differencing uses a fixed map when compared to running difference.
"""

import matplotlib.pyplot as plt
import matplotlib.colors as colors

import sunpy.map
from import sunpy.physics.differential_rotation import differential_rotate
from sunpy.data.sample import AIA_193_CUTOUT01_IMAGE, AIA_193_CUTOUT02_IMAGE, AIA_193_CUTOUT03_IMAGE

###############################################################################
# We create the MapSequence using the AIA_193_CUTOUT sample data.
# To create a MapSequence, we can call Map directly and add in a keyword to output a MapSequence instead.
aiamapseq = sunpy.map.Map(AIA_193_CUTOUT01_IMAGE, AIA_193_CUTOUT02_IMAGE,
                          AIA_193_CUTOUT03_IMAGE, sequence=True)

############################################################################
# In case of running difference, we loop through all the maps in the
# aiamapseq and differentially rotate each map in the MapSequence
# with respect to the previous map
# while in case of base difference we differentially
# rotate each map in the MapSequence to the time of the base map.
# We then store all the difference maps in a list.
base_diffmap = []
running_diffmap = []

# Storage for the differentially rotated maps. The first map is
# not differentially rotated.  We will also need the top left

# Let's do a base difference
drm = [aiamapseq[0]]
for m in aiamapseq[1:]:
    # differentially rotate each map and store it
    drm.append(differential_rotate(m, time=aiamapseq[0].date))

# Find the area common to all the differentially rotated maps
common_submap = drm[0]
for m in drm[1:]:
    bl, tr = overlap_coordinates(common_submap, m)
    common_submap = m.submap(bl, tr)
common_submap_dimensions = common_submap.dimensions


# Now do a base difference using the size of the common submap data.
base_difference_maps = [drm[0].submap(bl, tr)]
for m in drm[1:]:
    m_submap = m.submap(bl, tr)
    base_difference_data = base_difference_maps[0].data - m_submap.data
    base_difference_map = sunpy.map.Map(base_difference_data, m.meta)
    base_difference_map.plot_settings['cmap'] = plt.get_cmap('Greys_r')
    base_difference_map.plot_settings['norm'] = colors.LogNorm(100, base_difference_data.max())
    base_difference_maps.append(base_difference_map)


############################################################################
# This plots the original MapSequence
aiamapseq.peek()

##############################################################################
# This plots the final MapSequence obtained after implementing the base difference
base_difference_mapsequence = sunpy.map.MapSequence(base_difference_maps)
base_difference_mapsequence.peek()

############################################################################
# This plots the final MapSequence after implementing the running difference
result_mapseq = sunpy.map.MapSequence(running_diffmap)
result_mapseq.peek()
plt.show()
