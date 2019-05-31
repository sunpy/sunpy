# -*- coding: utf-8 -*-
"""
========================
Creating a Composite map
========================

How to create a composite map and use it to overplot two maps to
compare features.
"""
import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord
import astropy.units as u

import sunpy.map
import sunpy.data.sample

###############################################################################
# We start with the sample data. HMI shows the line-of-sight magnetic field at
# the photosphere while AIA 171 images show the resulting magnetic fields
# filled with hot plasma above, in the corona. We want to see what coronal
# features overlap with regions of strong line-of-sight magnetic fields.
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
hmi_map = sunpy.map.Map(sunpy.data.sample.HMI_LOS_IMAGE)

bottom_left = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=aia_map.coordinate_frame)
top_right = SkyCoord(800 * u.arcsec, 800 * u.arcsec, frame=aia_map.coordinate_frame)
aia_smap = aia_map.submap(bottom_left, top_right)
hmi_smap = hmi_map.submap(bottom_left, top_right)


##############################################################################
# Let's create a `~sunpy.map.CompositeMap` which includes both maps.
comp_map = sunpy.map.Map(aia_smap, hmi_smap, composite=True)

# Let's set the contours of the HMI map, the second image in our composite map
# (therefore the index is 1), from a few hundred to a thousand Gauss which
# is the typical field associated with umbral regions of active regions.
comp_map.set_levels(index=1, levels=[-1000, -500, -250, 250, 500, 1000])

##############################################################################
# Now let us look at the result. Notice that we can see the coronal structures
# present on the AIA image and how they correspond to the line of sight
# magnetic field.
fig = plt.figure()
ax = plt.subplot(projection=aia_smap)
comp_map.plot()
plt.show()
