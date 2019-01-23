# -*- coding: utf-8 -*-
"""
========================
Creating a Composite map
========================

In this example we make a composite map out of images provided by AIA and HMI.
AIA images look at the solar corona at a given temperature while an HMI image will
measure the line of sight magnetic field at the photospheric level.
"""

# Start by importing the necessary modules.
import sunpy.map
import sunpy.data.sample

##############################################################################
# Sunpy's sample data contains a selection of data
# which we will use for this example

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
hmi_map = sunpy.map.Map(sunpy.data.sample.HMI_LOS_IMAGE)

##############################################################################
# Now we will create the composite map which is a set of images stacked upon
# each other. To do this, we have to call `~sunpy.map.Map` with
# an extra keyword.

comp_map = sunpy.map.Map(aia_map, hmi_map, composite=True)

##############################################################################
# Now we want to see what coronal features overlap with regions of strong
# line of sight magnetic field. So we need to set the contour levels of
# our composite map using `.set_levels`.

# We want to contour the HMI map, which is the second image in our composite map.
# Therefore the index is 1.
# We will filter contours ranging from a few hundred to a thousand Gauss which
# is the typical field associated to umbral regions of Active Regions.
comp_map.set_levels(index=1, levels=[-1000, -500, -250, 250, 500, 1000])

##############################################################################
# Now let us look at the result. Notice that we can see the coronal structures
# present on the AIA image and how they correspond to the line of sight
# magnetic field.

comp_map.peek()
