"""
================================================================
Solar Orbiter Multi-Instrument Data Analysis and Coordination
================================================================

This example demonstrates advanced analysis techniques using Solar Orbiter data,
including combining EUI images with other Solar Orbiter instruments and
performing coordinate transformations for the non-Earth observer perspective.

Solar Orbiter's unique orbital perspective allows for novel observations of
solar phenomena from different viewpoints than Earth-based observatories.
"""
# sphinx_gallery_thumbnail_number = 3
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, AsinhStretch, LinearStretch
from astropy.time import Time

import sunpy.map
from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# First, let's demonstrate how to find Solar Orbiter's position relative to Earth
# and the Sun at a given time. This is important for understanding the
# observational context.

observation_time = Time('2022-03-30T10:00:00')

# Get Solar Orbiter's position
try:
    solo_coord = get_body_heliographic_stonyhurst('Solar Orbiter', observation_time)
    earth_coord = get_body_heliographic_stonyhurst('Earth', observation_time)

    print("Solar Orbiter Orbital Information")
    print("=" * 40)
    print(f"Observation Time: {observation_time}")
    print(f"Solar Orbiter Position: {solo_coord}")
    print(f"Earth Position: {earth_coord}")
    print(f"Distance from Sun: {solo_coord.radius:.3f}")
    print(f"Heliographic Longitude: {solo_coord.lon:.1f}")
    print(f"Heliographic Latitude: {solo_coord.lat:.1f}")

    # Calculate separation angle between Solar Orbiter and Earth as seen from Sun
    separation = solo_coord.separation(earth_coord)
    print(f"Angular separation from Earth: {separation:.1f}")

except Exception as e:
    print(f"Could not retrieve ephemeris data: {e}")
    print("Using nominal values for demonstration")
    solo_coord = SkyCoord(lon=45*u.deg, lat=5*u.deg, radius=0.8*u.AU,
                          frame='heliographic_stonyhurst',
                          obstime=observation_time)

###############################################################################
# Now let's search for and analyze different types of Solar Orbiter data.
# We'll look for both EUI images and magnetometer data from the same time period.

try:
    import sunpy_soar  # noqa: F401

    # Search for EUI data
    time_range = a.Time('2022-03-30T09:30:00', '2022-03-30T10:30:00')

    # EUI Full Sun Imager data
    eui_result = Fido.search(time_range,
                             a.Instrument('EUI'),
                             a.Level(2),
                             a.soar.ProductType('EUI-FSI174-IMAGE'))

    # MAG (magnetometer) data for context
    mag_result = Fido.search(time_range,
                             a.Instrument('MAG'),
                             a.Level(2))

    print(f"\nFound {len(eui_result)} EUI images")
    print(f"Found {len(mag_result)} MAG files")

    if len(eui_result) > 0:
        eui_files = Fido.fetch(eui_result[0])
        eui_map = sunpy.map.Map(eui_files[0])
        use_real_data = True
    else:
        use_real_data = False

except (ImportError, Exception):
    use_real_data = False
    print("Using simulated data for analysis demonstration")

if not use_real_data:
    # Create simulated EUI data based on AIA data
    from sunpy.data.sample import AIA_171_IMAGE
    aia_map = sunpy.map.Map(AIA_171_IMAGE)

    # Modify metadata to simulate Solar Orbiter EUI observation
    eui_meta = aia_map.meta.copy()
    eui_meta['instrume'] = 'EUI'
    eui_meta['obsrvtry'] = 'Solar Orbiter'
    eui_meta['detector'] = 'FSI'
    eui_meta['wavelnth'] = 174
    eui_meta['waveunit'] = 'Angstrom'

    # Simulate Solar Orbiter's perspective by adjusting observer coordinates
    eui_meta['hgln_obs'] = 45.0  # Heliographic longitude
    eui_meta['hglt_obs'] = 5.0   # Heliographic latitude
    eui_meta['dsun_obs'] = 0.8 * u.AU.to(u.m)  # Distance from Sun

    eui_map = sunpy.map.Map(aia_map.data, eui_meta)

###############################################################################
# Let's analyze the image to identify interesting solar features.
# We'll look for bright regions that might indicate active regions or flares.

# Calculate some basic image statistics
mean_intensity = np.mean(eui_map.data)
std_intensity = np.std(eui_map.data)
max_intensity = np.max(eui_map.data)

# Define threshold for "bright" regions (3 sigma above mean)
bright_threshold = mean_intensity + 3 * std_intensity

# Create a mask for bright regions
bright_mask = eui_map.data > bright_threshold

print(f"\nImage Analysis Results")
print("=" * 30)
print(f"Mean intensity: {mean_intensity:.1f}")
print(f"Standard deviation: {std_intensity:.1f}")
print(f"Maximum intensity: {max_intensity:.1f}")
print(f"Bright threshold (3σ): {bright_threshold:.1f}")
print(f"Percentage of bright pixels: {np.sum(bright_mask)/bright_mask.size*100:.2f}%")

###############################################################################
# Create a comprehensive analysis plot showing multiple views of the data

fig = plt.figure(figsize=(16, 12))

# Main EUI image
ax1 = fig.add_subplot(2, 3, 1, projection=eui_map)
norm1 = ImageNormalize(vmin=eui_map.data.min(), vmax=eui_map.data.max(),
                       stretch=AsinhStretch(0.01))
im1 = eui_map.plot(axes=ax1, norm=norm1, cmap='sdoaia171')
eui_map.draw_limb(axes=ax1, color='white', linewidth=2)
ax1.set_title(f'Solar Orbiter EUI {eui_map.wavelength}\nObserver at {getattr(eui_map, "hgln_obs", "unknown")}° longitude')

# Histogram of intensity values
ax2 = fig.add_subplot(2, 3, 2)
ax2.hist(eui_map.data.flatten(), bins=50, alpha=0.7, color='blue', density=True)
ax2.axvline(mean_intensity, color='red', linestyle='--', label='Mean')
ax2.axvline(bright_threshold, color='orange', linestyle='--', label='3σ threshold')
ax2.set_xlabel('Intensity')
ax2.set_ylabel('Density')
ax2.set_title('Intensity Distribution')
ax2.legend()
ax2.set_yscale('log')

# Bright regions overlay
ax3 = fig.add_subplot(2, 3, 3, projection=eui_map)
eui_map.plot(axes=ax3, norm=norm1, cmap='gray')
bright_overlay = np.ma.masked_where(~bright_mask, eui_map.data)
ax3.imshow(bright_overlay, origin='lower', extent=eui_map.wcs.array_index_to_world_values([0, 0], [eui_map.data.shape[1]-1, eui_map.data.shape[0]-1]).flatten(),
           cmap='Reds', alpha=0.7)
eui_map.draw_limb(axes=ax3, color='white', linewidth=2)
ax3.set_title('Bright Regions (>3σ)')

# Radial intensity profile
ax4 = fig.add_subplot(2, 3, 4)
# Calculate radial profile from disk center
y_center, x_center = np.array(eui_map.data.shape) / 2
y, x = np.ogrid[:eui_map.data.shape[0], :eui_map.data.shape[1]]
r = np.sqrt((x - x_center)**2 + (y - y_center)**2)

# Bin the data radially
r_max = min(x_center, y_center)
r_bins = np.linspace(0, r_max, 50)
radial_profile = []
for i in range(len(r_bins)-1):
    mask = (r >= r_bins[i]) & (r < r_bins[i+1])
    if np.any(mask):
        radial_profile.append(np.mean(eui_map.data[mask]))
    else:
        radial_profile.append(0)

r_centers = (r_bins[:-1] + r_bins[1:]) / 2
# Convert pixel radius to arcseconds (approximate)
pixel_scale = 4.4  # arcsec/pixel for EUI FSI (approximate)
r_arcsec = r_centers * pixel_scale

ax4.plot(r_arcsec, radial_profile, 'b-', linewidth=2)
ax4.set_xlabel('Distance from center [arcsec]')
ax4.set_ylabel('Mean intensity')
ax4.set_title('Radial Intensity Profile')
ax4.grid(True, alpha=0.3)

# Solar Orbiter orbit visualization
ax5 = fig.add_subplot(2, 3, 5)
# Create a simple orbital plot
orbit_angles = np.linspace(0, 2*np.pi, 100)
earth_orbit = np.ones_like(orbit_angles)  # 1 AU circle for Earth
solo_radius = getattr(solo_coord, 'radius', 0.8*u.AU).to(u.AU).value

# Plot orbits
ax5.plot(earth_orbit * np.cos(orbit_angles), earth_orbit * np.sin(orbit_angles),
         'b-', label='Earth orbit (1 AU)', linewidth=2)
ax5.plot(0, 0, 'yo', markersize=15, label='Sun')

# Plot current positions
earth_angle = 0  # Reference position
solo_angle = np.radians(getattr(solo_coord, 'lon', 45*u.deg).to(u.deg).value)
ax5.plot(1 * np.cos(earth_angle), 1 * np.sin(earth_angle),
         'bo', markersize=10, label='Earth')
ax5.plot(solo_radius * np.cos(solo_angle), solo_radius * np.sin(solo_angle),
         'ro', markersize=10, label='Solar Orbiter')

ax5.set_xlim(-1.5, 1.5)
ax5.set_ylim(-1.5, 1.5)
ax5.set_aspect('equal')
ax5.grid(True, alpha=0.3)
ax5.legend()
ax5.set_title('Orbital Positions')
ax5.set_xlabel('Distance [AU]')
ax5.set_ylabel('Distance [AU]')

# Feature detection summary
ax6 = fig.add_subplot(2, 3, 6)
ax6.axis('off')
feature_text = f"""
Solar Orbiter Analysis Summary
{'='*30}

Observation Details:
• Time: {eui_map.date.strftime('%Y-%m-%d %H:%M:%S')}
• Observer longitude: {getattr(eui_map, 'hgln_obs', 'Unknown')}°
• Observer latitude: {getattr(eui_map, 'hglt_obs', 'Unknown')}°
• Distance from Sun: {solo_radius:.2f} AU

Image Statistics:
• Dimensions: {eui_map.data.shape[0]}×{eui_map.data.shape[1]} pixels
• Mean intensity: {mean_intensity:.1f}
• Max intensity: {max_intensity:.1f}
• Bright pixels (>3σ): {np.sum(bright_mask)/bright_mask.size*100:.1f}%

Unique Perspective:
• Solar Orbiter's off-Earth-line view
• Enables stereoscopic observations
• Different viewing angle of solar features
"""
ax6.text(0.05, 0.95, feature_text, transform=ax6.transAxes, fontsize=10,
         verticalalignment='top', fontfamily='monospace')

plt.tight_layout()
plt.show()

###############################################################################
# Demonstrate coordinate transformation for off-limb features.
# Solar Orbiter's different perspective allows observation of features
# not visible from Earth.

print("\n" + "="*60)
print("Coordinate Transformation Analysis")
print("="*60)

# Define some example coordinates on the solar disk as seen from Solar Orbiter
example_coords = [
    SkyCoord(0*u.arcsec, 0*u.arcsec, frame=eui_map.coordinate_frame),  # Disk center
    SkyCoord(500*u.arcsec, 0*u.arcsec, frame=eui_map.coordinate_frame),  # East limb
    SkyCoord(-500*u.arcsec, 0*u.arcsec, frame=eui_map.coordinate_frame),  # West limb
    SkyCoord(0*u.arcsec, 500*u.arcsec, frame=eui_map.coordinate_frame),  # North limb
]

coord_names = ['Disk Center', 'East Limb', 'West Limb', 'North Limb']

print("Coordinates as seen from Solar Orbiter vs. Earth perspective:")
print("-" * 60)

for coord, name in zip(example_coords, coord_names):
    try:
        # Transform to heliographic coordinates
        hg_coord = coord.transform_to('heliographic_stonyhurst')
        print(f"{name:12}: Lat={hg_coord.lat.to(u.deg):7.1f}, Lon={hg_coord.lon.to(u.deg):7.1f}")
    except Exception:
        print(f"{name:12}: Coordinate transformation not available")

###############################################################################
# Finally, let's create a comparison plot showing how the same solar feature
# might look different from Earth vs. Solar Orbiter perspective.

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                               subplot_kw={'projection': eui_map})

# Plot 1: Solar Orbiter perspective (actual data)
norm = ImageNormalize(vmin=eui_map.data.min(), vmax=eui_map.data.max(),
                      stretch=AsinhStretch(0.01))
eui_map.plot(axes=ax1, norm=norm, cmap='sdoaia171')
eui_map.draw_limb(axes=ax1, color='white', linewidth=2)
eui_map.draw_grid(axes=ax1, color='white', alpha=0.3)

# Add annotation for unique viewing angle
ax1.text(0.02, 0.98, f'Solar Orbiter View\\nLongitude: {getattr(eui_map, "hgln_obs", "45")}°',
         transform=ax1.transAxes, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='black', alpha=0.7),
         color='white', fontweight='bold')

ax1.set_title('Solar Orbiter EUI Perspective')

# Plot 2: Simulated Earth perspective (rotated view)
# This is a simplified demonstration - real coordinate transformation would be more complex
rotated_data = np.fliplr(eui_map.data)  # Simple flip to simulate different viewing angle
earth_meta = eui_map.meta.copy()
earth_meta['hgln_obs'] = 0.0  # Earth's longitude
earth_meta['hglt_obs'] = 0.0  # Earth's latitude

earth_view_map = sunpy.map.Map(rotated_data, earth_meta)
earth_view_map.plot(axes=ax2, norm=norm, cmap='sdoaia171')
earth_view_map.draw_limb(axes=ax2, color='white', linewidth=2)
earth_view_map.draw_grid(axes=ax2, color='white', alpha=0.3)

ax2.text(0.02, 0.98, 'Earth-like View\\nLongitude: 0°',
         transform=ax2.transAxes, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='black', alpha=0.7),
         color='white', fontweight='bold')

ax2.set_title('Simulated Earth Perspective')

plt.tight_layout()
plt.show()

print("\nAnalysis Complete!")
print("This example demonstrates the unique scientific value of Solar Orbiter's")
print("orbital perspective for solar physics research.")
