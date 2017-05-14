# -*- coding: utf-8 -*-
"""
=================================
Drawing AIA Limb on STEREO Images
=================================


In this example we use a STEREO-B and an SDO image to demonstrate how to
overplot the limb as seen by AIA on an EUVI-B image. This makes use of
functionality added in Astropy 1.3.
"""

##############################################################################
# Start by importing the necessary modules.

from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import sunpy.coordinates
import sunpy.coordinates.wcs_utils
from sunpy.net import vso


##############################################################################
# The first step is to download some data, we are going to get an image from
# early 2011 when the STEREO spacecraft were roughly 90 deg seperated from the
# Earth.

stereo = (vso.attrs.Source('STEREO_B') &
          vso.attrs.Instrument('EUVI') &
          vso.attrs.Time('2011-01-01', '2011-01-01T00:10:00'))

aia = (vso.attrs.Instrument('AIA') &
       vso.attrs.Sample(24 * u.hour) &
       vso.attrs.Time('2011-01-01', '2011-01-02'))

wave = vso.attrs.Wavelength(30 * u.nm, 31 * u.nm)


vc = vso.VSOClient()
res = vc.query(wave, aia | stereo)


##############################################################################
# The results from VSO query:

print(res)


##############################################################################
# Download the files:

files = vc.get(res).wait()


##############################################################################
# Create a dictionary with the two maps, cropped down to full disk.

maps = {m.detector: m.submap((-1100, 1100) * u.arcsec,
                             (-1100, 1100) * u.arcsec) for m in sunpy.map.Map(files)}


##############################################################################
# Calculate points on the limb in the AIA image for the half that can be seen
# from STEREO.

r = maps['AIA'].rsun_obs - 1 * u.arcsec # remove the one arcsec so it's on disk.
# Adjust the following range if you only want to plot on STEREO_A
th = np.linspace(-180*u.deg, 0*u.deg)
x = r * np.sin(th)
y = r * np.cos(th)

coords = SkyCoord(x, y, frame=maps['AIA'].coordinate_frame)


##############################################################################
# Plot both maps

fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(1, 2, 1, projection=maps['AIA'])
maps['AIA'].plot(axes=ax1)
maps['AIA'].draw_limb()

ax2 = fig.add_subplot(1, 2, 2, projection=maps['EUVI'])
maps['EUVI'].plot(axes=ax2)
ax2.plot_coord(coords, color='w')
