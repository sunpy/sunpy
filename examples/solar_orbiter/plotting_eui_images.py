"""
========================================
Plotting Solar Orbiter EUI Images
========================================

This example demonstrates how to plot Solar Orbiter Extreme Ultraviolet Imager (EUI)
data with proper styling and coordinate information.

EUI has three telescopes: the Full Sun Imager (FSI) that images the whole Sun
in 174 Å and 304 Å, and two High Resolution Imagers (HRI) that image smaller
patches in 174 Å (EUV) and 1216 Å (Lyman-alpha).
"""
# sphinx_gallery_thumbnail_number = 2
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, AsinhStretch

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# For this example, we'll use sample data. In a real scenario, you would
# download EUI data using sunpy-soar as shown in the previous example.
#
# Let's create a mock Solar Orbiter EUI map to demonstrate plotting techniques.

try:
    import sunpy_soar  # noqa: F401

    # Search for actual EUI data
    time_range = a.Time('2022-03-30T10:00:00', '2022-03-30T10:30:00')
    instrument = a.Instrument('EUI')
    level = a.Level(2)
    product_type = a.soar.ProductType('EUI-FSI174-IMAGE')

    result = Fido.search(time_range, instrument, level, product_type)

    if len(result) > 0:
        files = Fido.fetch(result[0])
        eui_map = sunpy.map.Map(files[0])
        print(f"Using real EUI data: {eui_map.meta.get('filename', 'Unknown')}")
    else:
        raise ValueError("No data found")

except (ImportError, ValueError):
    print("Using simulated EUI data for demonstration")
    # For demonstration, let's use AIA data as a proxy and modify its metadata
    # to simulate EUI data characteristics
    import numpy as np
    from sunpy.data.sample import AIA_171_IMAGE

    # Load AIA map and modify it to simulate EUI characteristics
    aia_map = sunpy.map.Map(AIA_171_IMAGE)

    # Create simulated EUI metadata
    eui_meta = aia_map.meta.copy()
    eui_meta['instrume'] = 'EUI'
    eui_meta['obsrvtry'] = 'Solar Orbiter'
    eui_meta['detector'] = 'FSI'
    eui_meta['wavelnth'] = 174
    eui_meta['waveunit'] = 'Angstrom'
    eui_meta['rsun_ref'] = 696000000  # meters
    eui_meta['dsun_obs'] = 1.5e11  # meters (1 AU)

    # Create the simulated EUI map
    eui_map = sunpy.map.Map(aia_map.data, eui_meta)

###############################################################################
# Now let's create our first plot. EUI images benefit from asinh stretching
# to bring out both bright and faint features.

fig = plt.figure(figsize=(10, 6))

# First subplot: Basic plot with default settings
ax1 = fig.add_subplot(121, projection=eui_map)
eui_map.plot(axes=ax1)
eui_map.draw_limb(axes=ax1, color='white', linewidth=2)
ax1.set_title(f'EUI {eui_map.wavelength} - Default')

# Second subplot: Enhanced plot with asinh stretch
ax2 = fig.add_subplot(122, projection=eui_map)
norm = ImageNormalize(vmin=eui_map.data.min(), vmax=eui_map.data.max(),
                      stretch=AsinhStretch(0.01))
eui_map.plot(axes=ax2, norm=norm, cmap='sdoaia171')
eui_map.draw_limb(axes=ax2, color='white', linewidth=2)
eui_map.draw_grid(axes=ax2, color='white', alpha=0.5)
ax2.set_title(f'EUI {eui_map.wavelength} - Enhanced')

plt.tight_layout()
plt.show()

###############################################################################
# Let's create a more detailed plot showing coordinate information and
# adding some annotations.

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(projection=eui_map)

# Plot with enhanced normalization
norm = ImageNormalize(vmin=eui_map.data.min(), vmax=eui_map.data.max(),
                      stretch=AsinhStretch(0.01))
im = eui_map.plot(axes=ax, norm=norm, cmap='sdoaia171')

# Add solar limb and grid
eui_map.draw_limb(axes=ax, color='white', linewidth=2)
eui_map.draw_grid(axes=ax, color='white', alpha=0.3, linewidth=0.5)

# Add coordinate labels
ax.coords[0].set_axislabel('Solar X [arcsec]')
ax.coords[1].set_axislabel('Solar Y [arcsec]')

# Create a detailed title with observation information
title = (f"Solar Orbiter EUI {eui_map.detector} "
         f"{eui_map.wavelength} Image\\n"
         f"Observed: {eui_map.date.strftime('%Y-%m-%d %H:%M:%S')} UTC")
ax.set_title(title, fontsize=12, pad=20)

# Add a colorbar
cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Intensity [DN/s]', rotation=270, labelpad=20)

# Add some example annotations
# Mark the disk center
center = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=eui_map.coordinate_frame)
ax.plot_coord(center, 'x', color='red', markersize=10, markeredgewidth=2,
              markeredgecolor='white')
ax.text(0.02, 0.98, 'Disk Center', transform=ax.transAxes,
        verticalalignment='top', color='white', fontweight='bold')

plt.tight_layout()
plt.show()

###############################################################################
# Display some key information about the Solar Orbiter observation

print("\\n" + "="*50)
print("Solar Orbiter EUI Observation Details")
print("="*50)
print(f"Observatory: {eui_map.observatory}")
print(f"Instrument: {eui_map.instrument}")
print(f"Detector: {getattr(eui_map, 'detector', 'Unknown')}")
print(f"Wavelength: {eui_map.wavelength}")
print(f"Observation Time: {eui_map.date}")
print(f"Exposure Time: {eui_map.exposure_time}")
print(f"Image Dimensions: {eui_map.dimensions}")
print(f"Pixel Scale: {eui_map.scale}")
print(f"Observer Distance: {eui_map.dsun:.2e}")
print(f"Solar Radius: {eui_map.rsun_obs:.1f}")

###############################################################################
# EUI data often shows interesting coronal structures. Let's create a plot
# that highlights these features using different color schemes.

fig, axes = plt.subplots(1, 3, figsize=(18, 6),
                         subplot_kw={'projection': eui_map})

# Three different color schemes for EUI data
colormaps = ['sdoaia171', 'plasma', 'hot']
titles = ['SDO AIA Style', 'Plasma', 'Hot']

for ax, cmap, title in zip(axes, colormaps, titles):
    norm = ImageNormalize(vmin=eui_map.data.min(), vmax=eui_map.data.max(),
                          stretch=AsinhStretch(0.01))
    eui_map.plot(axes=ax, norm=norm, cmap=cmap)
    eui_map.draw_limb(axes=ax, color='white', linewidth=1)
    ax.set_title(f'EUI {eui_map.wavelength} - {title}')
    ax.coords[0].set_axislabel('')
    ax.coords[1].set_axislabel('')

# Only label the first subplot
axes[0].coords[1].set_axislabel('Solar Y [arcsec]')
axes[1].coords[0].set_axislabel('Solar X [arcsec]')

plt.tight_layout()
plt.show()