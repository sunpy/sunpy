# coding: utf-8
"""
================================================
Drawing and using a Great Arc and a Great Circle
================================================

How to define and draw a great arc on an image of the
Sun, and to extract intensity values along that arc.
"""
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.coordinates.utils import CoordinateVisibility, GreatArc
from sunpy.data.sample import AIA_171_IMAGE
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# We start with the sample data
m = sunpy.map.Map(AIA_171_IMAGE)

###############################################################################
# Let's define the start and end coordinates of the arc.
initial = SkyCoord(735 * u.arcsec, -471 * u.arcsec, frame=m.coordinate_frame)
target = SkyCoord(-100 * u.arcsec, 800 * u.arcsec, frame=m.coordinate_frame)

###############################################################################
# Create the great arc between the start and end points.
great_arc = GreatArc(initial, target)

###############################################################################
# Plot the great arc on the Sun.
fig = plt.figure()
ax = plt.subplot(projection=m)
m.plot(axes=ax)
plot_on_front = {"color": 'y', "linewidth": 5, "label": 'visible'}
plot_initial = {"color": 'r', "label": 'initial'}
plot_target = {"color": 'r', "label": 'target'}
ax.plot_coord(great_arc.coordinates, **plot_on_front)
ax.plot_coord(initial, 'x', **plot_initial)
ax.plot_coord(target, 'o', **plot_target)
plt.show()

###############################################################################
# Now we can calculate the nearest integer pixels of the data that correspond
# to the location of arc.
pixels = np.asarray(np.rint(m.world_to_pixel(great_arc.coordinates)), dtype=int)
x = pixels[0, :]
y = pixels[1, :]

###############################################################################
# Get the intensity along the arc from the start to the end point.
intensity_along_arc = m.data[y, x]

###############################################################################
# Define the angular location of each pixel along the arc from the start point
# to the end.
angles = great_arc.angles.to(u.deg)

###############################################################################
# Plot the intensity along the arc from the start to the end point.
fig, ax = plt.subplots()
ax.plot(angles, intensity_along_arc)
ax.set_xlabel('degrees of arc from start')
ax.set_ylabel('intensity')
ax.grid(linestyle='dotted')
plt.show()

###############################################################################
# The example above draws an arc along the inner angle directed from the start
# to the end coordinate.  The outer angle can also be used to define the arc.
great_arc = GreatArc(initial, target, use_inner_angle_direction=False)
ga = great_arc.coordinates
v = CoordinateVisibility(ga)

from_back_to_front = v.from_back_to_front
from_front_to_back = v.from_front_to_back
fig = plt.figure()
ax = plt.subplot(projection=m)
m.plot(axes=ax)
plot_on_back = {"color": 'c', "linewidth": 5, "label": 'not visible', "linestyle": ":"}
ax.plot_coord(ga[0:from_front_to_back], **plot_on_front)
ax.plot_coord(ga[from_front_to_back+1: from_back_to_front], **plot_on_back)
ax.plot_coord(ga[from_back_to_front+1:], **plot_on_front)
ax.plot_coord(initial, 'x', **plot_initial)
ax.plot_coord(target, 'o', **plot_target)
plt.show()


###############################################################################
# Great circles can also be drawn using the GreatArc object.  The following
# example creates a great circle that passes through two points on the solar
# surface, the first point seen from AIA and the second as seen from STEREO A.
stereo = (a.vso.Source('STEREO_B') &
          a.Instrument("EUVI") &
          a.Time('2011-01-01', '2011-01-01T00:10:00'))

aia = (a.Instrument('AIA') &
       a.Sample(24 * u.hour) &
       a.Time('2011-01-01', '2011-01-02'))

wave = a.Wavelength(30 * u.nm, 31 * u.nm)

res = Fido.search(wave, aia | stereo)
files = Fido.fetch(res)

if 'euvi' in files[0]:
    stereo_map = sunpy.map.Map(files[0])
    aia_map = sunpy.map.Map(files[1])
else:
    stereo_map = sunpy.map.Map(files[1])
    aia_map = sunpy.map.Map(files[0])

initial = SkyCoord(735 * u.arcsec, -471 * u.arcsec, frame=aia_map.coordinate_frame)
target = SkyCoord(-100 * u.arcsec, 800 * u.arcsec, frame=aia_map.coordinate_frame)

# Great circle as seen from AIA
great_circle = GreatArc(initial, target, points=1000, great_circle=True, use_inner_angle_direction=False)
c0 = great_circle.coordinates
c0v = CoordinateVisibility(c0)

# Great circle coordinates as seen from STEREO
c1 = c0.transform_to(stereo_map.coordinate_frame)
c1v = CoordinateVisibility(c1)

# The part of the arc which is visible from AIA and STEREO A.
both = np.logical_and(c0v.visible, c1v.visible)

# The part of the arc which is not visible from either AIA or STEREO A.
neither = np.logical_and(~c0v.visible, ~c1v.visible)

###############################################################################
# Determine the indices of the coordinates that are on and not on the
# observable disk of the Sun as seen from AIA and STEREO


def front_indices(visibility):
    """
    Helper function that analyzes a CoordinateVisibility object to return the
    indices of the coordinates that are visible from the observer coordinate.
    """
    n = len(visibility.visible)
    logic = np.roll(visibility.visible, -visibility.from_back_to_front - 1)
    indices = np.roll(np.arange(0, n), -visibility.from_back_to_front - 1)
    return indices[logic]


def back_indices(visibility):
    """
    Helper function that analyzes a CoordinateVisibility object to return the
    indices of the coordinates that are not visible from the observer coordinate.
    """
    n = len(visibility.visible)
    logic = np.roll(~visibility.visible, -visibility.from_front_to_back - 1)
    indices = np.roll(np.arange(0, n), -visibility.from_front_to_back - 1)
    return indices[logic]


c0_front_arc_indices = front_indices(c0v)
c0_back_arc_indices = back_indices(c0v)

c1_front_arc_indices = front_indices(c1v)
c1_back_arc_indices = back_indices(c1v)


###############################################################################
# Plot the great circle and its visibility on both the AIA and STEREO maps
fig = plt.figure(figsize=(10, 4))
ax1 = fig.add_subplot(1, 2, 1, projection=aia_map)
aia_map.plot(axes=ax1)

###############################################################################
# Set up the colors and linestyles we want and create the plot.
plot_visible_to_both = {"color": 'k', "label": 'visible to both'}
plot_visible_to_neither = {"color": 'r', "label": 'visible to neither'}

ax1.plot_coord(c0[c0_front_arc_indices], **plot_on_front)
ax1.plot_coord(c0[c0_back_arc_indices], **plot_on_back)
ax1.plot_coord(c0[both], **plot_visible_to_both)
ax1.plot_coord(c0[neither], **plot_visible_to_neither)
ax1.plot_coord(initial.transform_to(aia_map.coordinate_frame), 'x', **plot_initial)
ax1.plot_coord(target.transform_to(aia_map.coordinate_frame), 'o', **plot_target)
plt.legend()

ax2 = fig.add_subplot(1, 2, 2, projection=stereo_map)
stereo_map.plot(axes=ax2)
ax2.plot_coord(c1[c1_front_arc_indices], **plot_on_front)
ax2.plot_coord(c1[c1_back_arc_indices], **plot_on_back)
ax2.plot_coord(c1[both], **plot_visible_to_both)
ax2.plot_coord(c1[neither], **plot_visible_to_neither)
ax2.plot_coord(initial.transform_to(stereo_map.coordinate_frame), 'x', **plot_initial)
ax2.plot_coord(target.transform_to(stereo_map.coordinate_frame), 'o', **plot_target)

plt.show()
