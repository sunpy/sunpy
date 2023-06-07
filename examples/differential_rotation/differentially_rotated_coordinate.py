"""
====================================
Differentially rotating a coordinate
====================================

How to differentially rotate a coordinate.

The example uses the `~sunpy.coordinates.metaframes.RotatedSunFrame` coordinate
metaframe in `sunpy.coordinates` to apply differential rotation to a
coordinate.  See :ref:`sunpy-topic-guide-coordinates-rotatedsunframe` for more details on
using `~sunpy.coordinates.metaframes.RotatedSunFrame`.
"""
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.coordinates import RotatedSunFrame
from sunpy.data.sample import AIA_171_IMAGE

##############################################################################
# First, load an AIA observation and define a coordinate in its coordinate
# frame (here, helioprojective Cartesian).  The appropriate rate of rotation
# is determined from the heliographic latitude of the coordinate.

aiamap = sunpy.map.Map(AIA_171_IMAGE)
point = SkyCoord(187*u.arcsec, 283*u.arcsec, frame=aiamap.coordinate_frame)

##############################################################################
# We can differentially rotate this coordinate by using
# `~sunpy.coordinates.metaframes.RotatedSunFrame` with an array of observation
# times.  Let's define a daily cadence for +/- five days.

durations = np.concatenate([range(-5, 0), range(1, 6)]) * u.day
diffrot_point = SkyCoord(RotatedSunFrame(base=point, duration=durations))

##############################################################################
# To see what this coordinate looks like in "real" helioprojective
# Cartesian coordinates, we can transform it back to the original frame.
# Since these coordinates are represented in the original frame, they will not
# account for the changing position of the observer over this same time range.

transformed_diffrot_point = diffrot_point.transform_to(aiamap.coordinate_frame)
print(transformed_diffrot_point)

##############################################################################
# Let's plot the original coordinate and the differentially rotated
# coordinates on top of the AIA observation.

fig = plt.figure()
ax = fig.add_subplot(projection=aiamap)
aiamap.plot(axes=ax, clip_interval=(1., 99.95)*u.percent)

ax.plot_coord(point, 'ro', fillstyle='none', label='Original')
ax.plot_coord(transformed_diffrot_point, 'bo', fillstyle='none',
              label='Rotated')
ax.legend()

plt.show()
