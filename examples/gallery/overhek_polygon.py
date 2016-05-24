"""
=====================================================
Overplotting HEK feature/event polygons on SunPy maps
=====================================================

This example shows how to overplot HEK outlines on SunPy maps.
"""

##############################################################################
# Start by importing the necessary modules.

from __future__ import print_function, division
from datetime import timedelta
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import patches
from astropy import units as u

import sunpy.map
import sunpy.data.sample
from sunpy.net import hek
from sunpy.time import parse_time
from sunpy.physics.solar_rotation import rot_hpc

##############################################################################
# Load in an AIA map:
aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

##############################################################################
# Let's look for sunspots in the HEK close to the time of the AIA map. First
# create a client:
hek_client = hek.HEKClient()

##############################################################################
# Look for coronal holes detected using the SPoCA feature recognition method:
start_time = aia_map.date - timedelta(hours=2)
end_time = aia_map.date + timedelta(hours=2)
responses = hek_client.query(hek.attrs.Time(start_time, end_time), hek.attrs.CH, hek.attrs.FRM.Name == 'SPoCA')

##############################################################################
# Let's find the biggest coronal hole within 60 degrees north/south of the
# equator:
area = 0.0
for i, response in enumerate(responses):
    if response['area_atdiskcenter'] > area and np.abs(response['hgc_y']) < 60.0:
        area = response['area_atdiskcenter']
        response_index = i

##############################################################################
# Now let's get the boundary of the coronal hole
ch = responses[response_index]
p1 = ch["hpc_boundcc"][9: -2]
p2 = p1.split(',')
p3 = [v.split(" ") for v in p2]
ch_boundary = np.asarray([(eval(v[0]), eval(v[1])) for v in p3])
ch_date = parse_time(ch['event_starttime'])

##############################################################################
# The coronal hole was detected at a certain time.  To plot it on a map, we
# need to rotate it to the map observation time.
rotated_boundary = np.zeros_like(ch_boundary)
n = ch_boundary.shape[0]
for i in range(0, n):
    new_coords = rot_hpc(ch_boundary[i, 0] * u.arcsec,
                         ch_boundary[i, 1] * u.arcsec,
                         ch_date,
                         aia_map.date)
    rotated_boundary[i, 0] = new_coords[0].value
    rotated_boundary[i, 1] = new_coords[1].value

##############################################################################
# Now let's plot the rotated coronal hole boundary on the AIA map, and fill
# it with some matplotlib hatching.
fig = plt.figure()
ax = plt.subplot()
aia_map.plot()
coronal_hole = patches.Polygon(rotated_boundary, color='white', fill=False)
ax.set_title('{:s}\n{:s}'.format(aia_map.name, ch['frm_specificid']))
ax.add_artist(coronal_hole)
plt.colorbar()
plt.tight_layout()
