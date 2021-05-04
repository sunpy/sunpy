"""
===========================================
Drawing the AIA limb on a STEREO EUVI image
===========================================

In this example we use a STEREO-B and an SDO image to demonstrate how to
overplot the limb as seen by AIA on an EUVI-B image. Then we overplot the AIA
coordinate grid on the STEREO image.
"""
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.coordinates.wcs_utils
import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

##############################################################################
# The first step is to download some data, we are going to get an image from
# early 2011 when the STEREO spacecraft were roughly 90 deg separated from the
# Earth.

stereo = (a.Source('STEREO_B') &
          a.Instrument("EUVI") &
          a.Time('2011-01-01', '2011-01-01T00:10:00'))

aia = (a.Instrument.aia &
       a.Sample(24 * u.hour) &
       a.Time('2011-01-01', '2011-01-02'))

wave = a.Wavelength(30 * u.nm, 31 * u.nm)
result = Fido.search(wave, aia | stereo)

###############################################################################
# Let's inspect the result and download the files.

print(result)
downloaded_files = Fido.fetch(result)
print(downloaded_files)

##############################################################################
# Let's create a dictionary with the two maps, which we crop to full disk.

maps = {m.detector: m.submap(SkyCoord([-1100, 1100], [-1100, 1100],
                                      unit=u.arcsec, frame=m.coordinate_frame))
        for m in sunpy.map.Map(downloaded_files)}

##############################################################################
# Now, let's plot both maps, and we draw the limb as seen by AIA onto the
# EUVI image. We remove the part of the limb that is hidden because it is on
# the far side of the Sun from STEREO's point of view.

fig = plt.figure(figsize=(10, 4))
ax1 = fig.add_subplot(1, 2, 1, projection=maps['AIA'])
maps['AIA'].plot(axes=ax1)
maps['AIA'].draw_limb()

ax2 = fig.add_subplot(1, 2, 2, projection=maps['EUVI'])
maps['EUVI'].plot(axes=ax2)
visible, hidden = maps['AIA'].draw_limb()
hidden.remove()

##############################################################################
# Let's also plot the helioprojective coordinate grid as seen by SDO on the
# STEREO image.

fig = plt.figure()
ax = plt.subplot(projection=maps['EUVI'])

maps['EUVI'].plot()

# Move the title so it does not clash with the extra labels.
tx, ty = ax.title.get_position()
ax.title.set_position([tx, ty + 0.08])

# Change the default grid labels.
stereo_x, stereo_y = ax.coords
stereo_x.set_axislabel("Helioprojective Longitude (STEREO B) [arcsec]")
stereo_y.set_axislabel("Helioprojective Latitude (STEREO B) [arcsec]")

# Add a new coordinate overlay in the SDO frame.
overlay = ax.get_coords_overlay(maps['AIA'].coordinate_frame)
overlay.grid()

# Configure the grid:
x, y = overlay

# Set the ticks to be on the top and left axes.
x.set_ticks_position('tr')
y.set_ticks_position('tr')

# Wrap the longitude at 180 deg rather than the default 360.
x.set_coord_type('longitude', 180.)

# Change the defaults to arcseconds
x.set_major_formatter('s.s')
y.set_major_formatter('s.s')

# Add axes labels
x.set_axislabel("Helioprojective Longitude (SDO) [arcsec]")
y.set_axislabel("Helioprojective Latitude (SDO) [arcsec]")

plt.show()
