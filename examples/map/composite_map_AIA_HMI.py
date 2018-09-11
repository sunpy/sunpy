# -*- coding: utf-8 -*-
"""
=============
Composite map
=============
In this example we make a composite map out of images provided by AIA and HMI.
"""

##############################################################################
# Start by importing the necessary modules.

import sunpy.map
import sunpy.data.sample

##############################################################################
# Sunpy sample data contains a number of suitable maps, for this example we
# will use both AIA 171 and HMI magnetogram
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

hmi_map = sunpy.map.Map(sunpy.data.sample.HMI_LOS_IMAGE)

#creating the composite map object
comp_map = sunpy.map.Map(aia_map, hmi_map, composite = True)

#drawing the countours over the hmi_map note that the hmi_map have index = 1
#according to our definition of comp_map
comp_map.set_levels(index = 1, levels = [-1000,-500,-250,250,500,1000])

#having a look at the map
comp_map.peek()

#It is also possible to set up the alpha value for a layer in the CompositeMap.
comp_map.set_alpha(index = 1, alpha = 0.5)

#having another look at the map
comp_map.peek()

##############################################################################