# -*- coding: utf-8 -*-
"""
=============
Composite map
=============
In this example we make a composite map out of images provided by AIA and HMI.
AIA images are looking at the Sun's Corona intensity at a given temperature
while an HMI image will measure the magnetic field (in this case the 
Line-Of-Sight component) at photospheric level.

This procedure consist into taking a pair of images and plot features of one
of the images on the top of the other. Therefore, we can see on the coronal
structures present on AIA images what is the field associated to them. 

To filter what magnetic structures we want to see represented on our map we will
use the levels to define what contours will be draw which, in this case, it 
will be selected contours ranging from a few hundred to a thousand Gauss which
is the tipicall field associated to the umbral regions of Active Regions. If 
the contour values include numbers of the order 1e1 - 1e2 Gauss quiet Sun 
regions would also be ploted on the top of the original map which in practice 
would fill the original AIA map with magnetic field readings. 

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
comp_map = sunpy.map.Map(aia_map, hmi_map, composite=True)

#drawing the countours over the hmi_map note that the hmi_map have index = 1
#according to our definition of comp_map
comp_map.set_levels(index=1, levels=[-1000,-500,-250,250,500,1000])

#having a look at the map
comp_map.peek()

##############################################################################
