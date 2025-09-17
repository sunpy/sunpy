"""
=======================================================
Solar Orbiter Stereoscopic Solar Observations
=======================================================

This example demonstrates how to use Solar Orbiter's unique orbital perspective
for stereoscopic observations of solar features in coordination with Earth-based
observatories like SDO/AIA.

Stereoscopic observations enable 3D reconstruction of solar structures and
provide insights into the spatial distribution of coronal plasma that are
impossible to obtain from a single viewpoint.
"""
# sphinx_gallery_thumbnail_number = 2
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.visualization import ImageNormalize, AsinhStretch

import sunpy.map
from sunpy.coordinates import get_body_heliographic_stonyhurst
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# Define observation time and search for coordinated observations
# We'll look for observations that are close in time from both Solar Orbiter
# and Earth-based instruments.

observation_time = Time('2022-03-30T10:00:00')
time_window = 30 * u.minute  # Allow 30-minute window for coordination

print(f"Searching for stereoscopic observations around {observation_time}")

# Define search parameters
start_time = observation_time - time_window
end_time = observation_time + time_window
time_range = a.Time(start_time, end_time)

###############################################################################
# Calculate the orbital positions of Solar Orbiter and Earth

try:
    solo_coord = get_body_heliographic_stonyhurst('Solar Orbiter', observation_time)
    earth_coord = get_body_heliographic_stonyhurst('Earth', observation_time)

    separation_angle = solo_coord.separation(earth_coord)

    print(f"\\nOrbital Configuration:")
    print(f"Solar Orbiter: {solo_coord}")
    print(f"Earth: {earth_coord}")
    print(f"Separation angle: {separation_angle:.1f}")

    # Calculate baseline length for stereoscopy
    baseline_au = separation_angle.to(u.rad).value * solo_coord.radius.to(u.AU).value
    print(f"Stereoscopic baseline: {baseline_au:.3f} AU")

except Exception as e:
    print(f"Using nominal orbital positions: {e}")
    # Use nominal values for demonstration
    solo_coord = SkyCoord(lon=45*u.deg, lat=5*u.deg, radius=0.8*u.AU,
                          frame='heliographic_stonyhurst', obstime=observation_time)
    earth_coord = SkyCoord(lon=0*u.deg, lat=0*u.deg, radius=1.0*u.AU,
                           frame='heliographic_stonyhurst', obstime=observation_time)
    separation_angle = 45*u.deg

###############################################################################
# Search for coordinated observations

try:
    import sunpy_soar  # noqa: F401

    # Search for Solar Orbiter EUI data
    solo_query = Fido.search(time_range,
                             a.Instrument('EUI'),
                             a.Level(2),
                             a.soar.ProductType('EUI-FSI174-IMAGE'))

    print(f"Found {len(solo_query)} Solar Orbiter EUI images")

    # Search for SDO/AIA data from Earth perspective
    aia_query = Fido.search(time_range,
                            a.Instrument('AIA'),
                            a.Wavelength(171*u.angstrom),
                            a.Sample(10*u.minute))

    print(f"Found {len(aia_query)} SDO/AIA images")

    # Download data if available
    if len(solo_query) > 0 and len(aia_query) > 0:
        solo_files = Fido.fetch(solo_query[0])
        aia_files = Fido.fetch(aia_query[0])

        solo_map = sunpy.map.Map(solo_files[0])
        aia_map = sunpy.map.Map(aia_files[0])
        use_real_data = True
    else:
        use_real_data = False

except ImportError:
    print("sunpy-soar not available, using simulated data")
    use_real_data = False

###############################################################################
# Create simulated stereoscopic data if real data unavailable

if not use_real_data:
    from sunpy.data.sample import AIA_171_IMAGE

    # Load AIA data as baseline
    aia_map_original = sunpy.map.Map(AIA_171_IMAGE)

    # Create simulated Solar Orbiter view
    solo_meta = aia_map_original.meta.copy()
    solo_meta['instrume'] = 'EUI'
    solo_meta['obsrvtry'] = 'Solar Orbiter'
    solo_meta['detector'] = 'FSI'
    solo_meta['hgln_obs'] = 45.0  # 45 degrees longitude from Earth
    solo_meta['hglt_obs'] = 5.0   # 5 degrees latitude
    solo_meta['dsun_obs'] = 0.8 * u.AU.to(u.m)

    # Simulate perspective shift by rotating the image slightly
    # This is a simplified representation of the actual coordinate transformation
    rotated_data = np.roll(aia_map_original.data, shift=50, axis=1)
    solo_map = sunpy.map.Map(rotated_data, solo_meta)

    # Keep AIA as Earth perspective
    aia_map = aia_map_original

    print("Using simulated stereoscopic observations")

###############################################################################
# Create comparison plot showing both perspectives

fig = plt.figure(figsize=(16, 10))

# Plot 1: Solar Orbiter perspective
ax1 = fig.add_subplot(2, 3, 1, projection=solo_map)
norm1 = ImageNormalize(vmin=solo_map.data.min(), vmax=solo_map.data.max(),
                       stretch=AsinhStretch(0.01))
solo_map.plot(axes=ax1, norm=norm1, cmap='sdoaia171')
solo_map.draw_limb(axes=ax1, color='white', linewidth=2)
ax1.set_title(f'Solar Orbiter EUI 174 Å\\nLon: {getattr(solo_map, "hgln_obs", 45):.0f}°')

# Plot 2: Earth perspective (SDO/AIA)
ax2 = fig.add_subplot(2, 3, 2, projection=aia_map)
norm2 = ImageNormalize(vmin=aia_map.data.min(), vmax=aia_map.data.max(),
                       stretch=AsinhStretch(0.01))
aia_map.plot(axes=ax2, norm=norm2, cmap='sdoaia171')
aia_map.draw_limb(axes=ax2, color='white', linewidth=2)
ax2.set_title(f'SDO/AIA 171 Å\\nLon: 0°')

# Plot 3: Difference image to highlight perspective differences
ax3 = fig.add_subplot(2, 3, 3)

# Ensure both images have the same dimensions for comparison
if solo_map.data.shape == aia_map.data.shape:
    # Normalize both images to the same scale for comparison
    solo_norm = (solo_map.data - np.mean(solo_map.data)) / np.std(solo_map.data)
    aia_norm = (aia_map.data - np.mean(aia_map.data)) / np.std(aia_map.data)

    diff_image = solo_norm - aia_norm

    im3 = ax3.imshow(diff_image, origin='lower', cmap='RdBu_r',
                     vmin=-3, vmax=3)
    ax3.set_title('Perspective Difference\\n(Solar Orbiter - Earth)')
    plt.colorbar(im3, ax=ax3, label='Normalized Difference')
else:
    ax3.text(0.5, 0.5, 'Images have different\\ndimensions',
             transform=ax3.transAxes, ha='center', va='center')
    ax3.set_title('Perspective Difference')

ax3.set_xlabel('X [pixels]')
ax3.set_ylabel('Y [pixels]')

# Plot 4: Orbital configuration
ax4 = fig.add_subplot(2, 3, 4)
ax4.set_aspect('equal')

# Plot orbits
orbit_angles = np.linspace(0, 2*np.pi, 100)
earth_orbit = np.ones_like(orbit_angles)
ax4.plot(earth_orbit * np.cos(orbit_angles), earth_orbit * np.sin(orbit_angles),
         'b--', alpha=0.5, label='Earth orbit')

# Plot Sun
ax4.plot(0, 0, 'yo', markersize=15, label='Sun')

# Plot current positions
earth_x, earth_y = 1, 0  # Earth at reference position
solo_x = 0.8 * np.cos(np.radians(45))  # Solar Orbiter position
solo_y = 0.8 * np.sin(np.radians(45))

ax4.plot(earth_x, earth_y, 'bo', markersize=10, label='Earth')
ax4.plot(solo_x, solo_y, 'ro', markersize=10, label='Solar Orbiter')

# Draw lines of sight to show viewing angles
los_length = 1.5
ax4.arrow(earth_x, earth_y, -los_length, 0, head_width=0.05,
          head_length=0.05, fc='blue', ec='blue', alpha=0.7)
ax4.arrow(solo_x, solo_y, -los_length*np.cos(np.radians(45)),
          -los_length*np.sin(np.radians(45)), head_width=0.05,
          head_length=0.05, fc='red', ec='red', alpha=0.7)

ax4.set_xlim(-2, 1.5)
ax4.set_ylim(-1, 1.5)
ax4.grid(True, alpha=0.3)
ax4.legend()
ax4.set_title(f'Orbital Configuration\\nSeparation: {separation_angle:.0f}')
ax4.set_xlabel('Distance [AU]')
ax4.set_ylabel('Distance [AU]')

# Plot 5: Stereoscopic analysis - height estimation
ax5 = fig.add_subplot(2, 3, 5)

# Simulate height measurement using stereoscopic parallax
# This is a simplified demonstration of the principle
heights_rs = np.linspace(1.0, 2.0, 100)  # Heights in solar radii
baseline_rad = separation_angle.to(u.rad).value
distance_au = solo_coord.radius.to(u.AU).value

# Calculate apparent displacement due to parallax
# parallax = baseline * height / distance
parallax_arcsec = []
for h in heights_rs:
    height_au = h * 696000 / u.AU.to(u.km)  # Convert Rs to AU
    parallax_rad = baseline_rad * height_au / distance_au
    parallax_arcsec.append(parallax_rad * u.rad.to(u.arcsec))

ax5.plot(heights_rs, parallax_arcsec, 'g-', linewidth=2)
ax5.set_xlabel('Feature Height [R☉]')
ax5.set_ylabel('Parallax Displacement [arcsec]')
ax5.set_title('Stereoscopic Height Sensitivity')
ax5.grid(True, alpha=0.3)

# Highlight detectable range
detectable_threshold = 1.0  # arcsec (typical resolution limit)
ax5.axhline(detectable_threshold, color='red', linestyle='--',
            label=f'Detection limit ({detectable_threshold}" resolution)')
ax5.legend()

# Plot 6: Time series showing observation coordination
ax6 = fig.add_subplot(2, 3, 6)

# Create time series showing observation windows
time_hours = np.linspace(-2, 2, 100)  # Hours around target time
earth_visibility = np.ones_like(time_hours)  # Earth always sees near side
solo_visibility = np.ones_like(time_hours)   # Solar Orbiter sees from angle

# Add some realistic gaps (e.g., data downlink, instrument cycles)
earth_gaps = np.abs(time_hours) < 0.1  # Brief gap around target time
solo_gaps = (np.abs(time_hours - 0.5) < 0.2) | (np.abs(time_hours + 0.8) < 0.1)

earth_visibility[earth_gaps] = 0
solo_visibility[solo_gaps] = 0

ax6.fill_between(time_hours, 0, earth_visibility, alpha=0.5, color='blue',
                 label='Earth observations')
ax6.fill_between(time_hours, 1.1, 1.1 + solo_visibility, alpha=0.5, color='red',
                 label='Solar Orbiter observations')

# Mark coordination window
ax6.axvspan(-0.5, 0.5, alpha=0.2, color='green', label='Coordination window')
ax6.axvline(0, color='black', linestyle=':', label='Target time')

ax6.set_xlim(-2, 2)
ax6.set_ylim(-0.1, 2.3)
ax6.set_xlabel('Time [hours]')
ax6.set_ylabel('Observation availability')
ax6.set_title('Observation Coordination')
ax6.legend()
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

###############################################################################
# Analyze stereoscopic capabilities and science potential

print("\\n" + "="*60)
print("Stereoscopic Analysis Results")
print("="*60)

# Calculate key stereoscopic parameters
baseline_km = separation_angle.to(u.rad).value * solo_coord.radius.to(u.km).value
resolution_limit = 1.0 * u.arcsec  # Typical resolution

# Height sensitivity calculation
min_detectable_height = (resolution_limit.to(u.rad).value *
                        solo_coord.radius.to(u.km).value / baseline_km *
                        u.AU.to(u.km) / 696000)  # In solar radii

print(f"Stereoscopic Configuration:")
print(f"  Separation angle: {separation_angle:.1f}")
print(f"  Baseline length: {baseline_km:.0f} km")
print(f"  Solar Orbiter distance: {solo_coord.radius.to(u.AU):.2f}")

print(f"\\nHeight Sensitivity:")
print(f"  Minimum detectable height: {min_detectable_height:.2f} R☉")
print(f"  (assuming {resolution_limit} resolution limit)")

print(f"\\nScience Applications:")
print(f"  • 3D reconstruction of coronal loops")
print(f"  • Height determination of CME fronts")
print(f"  • Depth resolution of active region structures")
print(f"  • Validation of coronal magnetic field models")

# Identify optimal observation conditions
if separation_angle > 20*u.deg and separation_angle < 120*u.deg:
    print(f"  ✓ Excellent separation angle for stereoscopy")
elif separation_angle > 10*u.deg:
    print(f"  ✓ Good separation angle for stereoscopy")
else:
    print(f"  ⚠ Limited separation angle - reduced stereoscopic sensitivity")

if solo_coord.radius < 0.7*u.AU:
    print(f"  ✓ Close Solar Orbiter distance enhances resolution")

###############################################################################
# Create a summary plot showing the science potential

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left: 3D visualization concept
ax1.set_aspect('equal')

# Draw Sun
sun_circle = plt.Circle((0, 0), 0.1, color='yellow', label='Sun')
ax1.add_patch(sun_circle)

# Draw a simulated coronal loop in 3D (projected to 2D)
loop_angles = np.linspace(0, np.pi, 50)
loop_height = 0.3
loop_x = loop_height * np.cos(loop_angles)
loop_y = 0.1 + 0.05 * np.sin(2*loop_angles)  # Some structure
loop_z = loop_height * np.sin(loop_angles)

# Earth perspective (view from x-axis)
earth_x_view = loop_x
earth_y_view = loop_y
ax1.plot(earth_x_view, earth_y_view, 'b-', linewidth=3, label='Earth view')

# Solar Orbiter perspective (view from 45° angle)
angle_rad = np.radians(45)
solo_x_view = loop_x * np.cos(angle_rad) - loop_z * np.sin(angle_rad)
solo_y_view = loop_y
ax1.plot(solo_x_view, solo_y_view, 'r-', linewidth=3, label='Solar Orbiter view')

# Observer positions
ax1.plot(1.0, 0, 'bo', markersize=10, label='Earth')
ax1.plot(0.7*np.cos(angle_rad), 0.7*np.sin(angle_rad), 'ro', markersize=10,
         label='Solar Orbiter')

ax1.set_xlim(-0.6, 1.2)
ax1.set_ylim(-0.3, 0.6)
ax1.legend()
ax1.set_title('3D Coronal Structure Reconstruction')
ax1.set_xlabel('Solar X [R☉]')
ax1.set_ylabel('Solar Y [R☉]')

# Right: Height determination accuracy
ax2.set_aspect('equal')

separations = np.linspace(5, 180, 100)
height_accuracies = []

for sep in separations:
    sep_rad = np.radians(sep)
    baseline = sep_rad * 0.8  # AU (approximate Solar Orbiter distance)
    # Height accuracy inversely proportional to baseline
    accuracy = 1.0 / (baseline * 10)  # Simplified relationship
    height_accuracies.append(accuracy)

ax2.plot(separations, height_accuracies, 'g-', linewidth=3)
ax2.axvline(separation_angle.to(u.deg).value, color='red', linestyle='--',
            linewidth=2, label=f'Current separation ({separation_angle:.0f})')

ax2.set_xlabel('Separation Angle [degrees]')
ax2.set_ylabel('Relative Height Uncertainty')
ax2.set_title('Stereoscopic Height Accuracy')
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.set_xlim(0, 180)

plt.tight_layout()
plt.show()

print("\\nThis example demonstrates the unique capabilities of Solar Orbiter")
print("for stereoscopic solar physics, enabling 3D studies impossible from")
print("a single vantage point.")
