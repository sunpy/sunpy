"""
=========================================
Tracking and Co-aligning an Active Region
=========================================

This example demonstrates how to track a particular region as a function of time,
create a cutout around that region, and align it at each time step to get an aligned datacube.
"""

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.coordinates import propagate_with_solar_surface
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# `sunpy.net.Fido` is the primary interface for searching and downloading data.
# The following code performs a search for solar data using the Fido client from the sunpy.net module.
query = Fido.search(a.Time('2012/3/4', '2012/3/5'),
                    a.Instrument.aia,
                    a.Wavelength(171*u.angstrom),
                    a.Sample(360*u.minute))
print(query)

###############################################################################
# Let's inspect the results and download the files.
files = Fido.fetch(query)
print(files)

###############################################################################
# Create a sequence of maps from the downloaded files.

m_seq = sunpy.map.Map(files, sequence=True)

# Plot the sequence of maps.


fig = plt.figure(figsize=(24, 8))
for i, m in enumerate(m_seq):
    ax = fig.add_subplot(1, len(m_seq), i+1, projection=m)
    m.plot(axes=ax)
    ax.set_xlabel(' ')
    ax.set_ylabel(' ')
plt.subplots_adjust(wspace=0.3)


# The above images are the result of the query: 4 images at 6-hour intervals.

###############################################################################
# Define the region of interest using a SkyCoord.


corner = SkyCoord(Tx=-375*u.arcsec, Ty=0*u.arcsec, frame=m_seq[0].coordinate_frame)
print(m_seq[0].world_to_pixel(corner))

# Create a cutout around the defined region.
m_cutout = m_seq[0].submap(corner, width=500*u.arcsec, height=500*u.arcsec)

# Plot the cutout region.
fig = plt.figure()
ax = fig.add_subplot(projection=m_cutout)
m_cutout.plot(axes=ax)
plt.show()

###############################################################################
# Track and co-align the region across the sequence of maps using solar rotation.


fig = plt.figure(figsize=(24, 8))
for i, m in enumerate(m_seq):
    ax = fig.add_subplot(1, len(m_seq), i+1, projection=m)
    m.plot(axes=ax)
    ax.set_xlabel(' ')
    ax.set_ylabel(' ')
    plt.subplots_adjust(wspace=0.3)
    with propagate_with_solar_surface():
        blc = m_cutout.bottom_left_coord.transform_to(m.coordinate_frame)
        trc = m_cutout.top_right_coord.transform_to(m.coordinate_frame)
        m.draw_quadrangle(blc, top_right=trc)

###############################################################################
# Aligns all images to the World Coordinate System (WCS) of the cutout.
# `propagate_with_solar_surface` is a context manager that adjusts for solar rotation,
# ensuring that regions on the Sun's surface remain in the correct position as the Sun rotates.

with propagate_with_solar_surface():
    m_seq_aligned = sunpy.map.Map([m.reproject_to(m_cutout.wcs) for m in m_seq], sequence=True)

# `reproject_to` reprojects each map in the sequence to the WCS of the cutout map,
# aligning all images to the same reference frame.

# Plot the aligned sequence of maps.
fig = plt.figure(figsize=(24, 8))
for i, m in enumerate(m_seq_aligned):
    ax = fig.add_subplot(1, len(m_seq_aligned), i+1, projection=m)
    m.plot(axes=ax, cmap='sdoaia171', title=m_seq[i].date)
    ax.set_xlabel(' ')
    ax.set_ylabel(' ')
plt.subplots_adjust(wspace=0.3)
# This code aligns a sequence of solar images to a common reference frame, taking into account solar rotation,
# and then plots them side by side in a single figure. Each subplot shows one map from the aligned sequence.
