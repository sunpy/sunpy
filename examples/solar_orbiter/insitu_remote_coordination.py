"""
=======================================================
Solar Orbiter In-Situ and Remote Sensing Coordination
=======================================================

This example demonstrates how to coordinate Solar Orbiter's in-situ measurements
(like magnetic field data) with remote sensing observations (like EUI images).
This coordination is one of Solar Orbiter's unique scientific capabilities.

Solar Orbiter carries both remote sensing instruments that observe the Sun from
a distance and in-situ instruments that measure the local solar wind properties
at the spacecraft location.
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter

import astropy.units as u
from astropy.time import Time
from astropy.visualization import ImageNormalize, AsinhStretch

import sunpy.map
import sunpy.timeseries
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# Define the time range for our analysis. We'll look for both remote sensing
# and in-situ data from the same period.

start_time = '2022-03-30T06:00:00'
end_time = '2022-03-30T18:00:00'
time_range = a.Time(start_time, end_time)

print(f"Searching for Solar Orbiter data from {start_time} to {end_time}")

###############################################################################
# Search for different types of Solar Orbiter data using sunpy-soar

try:
    import sunpy_soar  # noqa: F401

    # Search for EUI images (remote sensing)
    eui_query = Fido.search(time_range,
                            a.Instrument('EUI'),
                            a.Level(2),
                            a.soar.ProductType('EUI-FSI174-IMAGE'))

    # Search for MAG data (in-situ magnetic field)
    mag_query = Fido.search(time_range,
                            a.Instrument('MAG'),
                            a.Level(2))

    # Search for SWA data (in-situ plasma measurements)
    swa_query = Fido.search(time_range,
                            a.Instrument('SWA'),
                            a.Level(2))

    print(f"Found {len(eui_query)} EUI images")
    print(f"Found {len(mag_query)} MAG files")
    print(f"Found {len(swa_query)} SWA files")

    # Download a subset of the data
    if len(eui_query) > 0:
        eui_files = Fido.fetch(eui_query[0:2])  # Download first 2 EUI files
    if len(mag_query) > 0:
        mag_files = Fido.fetch(mag_query[0:1])   # Download first MAG file

    data_available = len(eui_query) > 0 or len(mag_query) > 0

except ImportError:
    print("sunpy-soar not available. Using simulated data for demonstration.")
    data_available = False

###############################################################################
# If real data isn't available, create simulated datasets for demonstration

if not data_available:
    # Create simulated EUI data
    from sunpy.data.sample import AIA_171_IMAGE
    aia_map = sunpy.map.Map(AIA_171_IMAGE)

    # Modify to simulate Solar Orbiter EUI
    eui_meta = aia_map.meta.copy()
    eui_meta['instrume'] = 'EUI'
    eui_meta['obsrvtry'] = 'Solar Orbiter'
    eui_meta['detector'] = 'FSI'
    eui_meta['hgln_obs'] = 45.0  # Solar Orbiter longitude
    eui_meta['hglt_obs'] = 5.0   # Solar Orbiter latitude
    eui_meta['dsun_obs'] = 0.8 * u.AU.to(u.m)

    eui_map = sunpy.map.Map(aia_map.data, eui_meta)

    # Create simulated magnetometer data
    time_obs = Time(start_time) + np.linspace(0, 12, 720) * u.hour

    # Simulate magnetic field components (RTN coordinates)
    # RTN = Radial, Tangential, Normal coordinates
    np.random.seed(42)  # For reproducible results

    # Typical solar wind magnetic field values
    b_r = 3.0 + 2.0 * np.sin(2*np.pi*np.arange(720)/720) + 0.5 * np.random.randn(720)
    b_t = 1.0 + 1.5 * np.cos(2*np.pi*np.arange(720)/720) + 0.3 * np.random.randn(720)
    b_n = 0.5 * np.sin(4*np.pi*np.arange(720)/720) + 0.2 * np.random.randn(720)
    b_total = np.sqrt(b_r**2 + b_t**2 + b_n**2)

    # Create synthetic timeseries data
    mag_data = {
        'B_R': b_r,
        'B_T': b_t,
        'B_N': b_n,
        'B_TOTAL': b_total
    }

    # Simulate plasma data
    proton_density = 5.0 + 3.0 * np.exp(-((np.arange(720)-360)/100)**2) + 0.5 * np.random.randn(720)
    proton_speed = 400 + 100 * np.sin(2*np.pi*np.arange(720)/720) + 10 * np.random.randn(720)
    proton_temp = 50000 + 20000 * np.random.randn(720)

    plasma_data = {
        'N_P': np.maximum(proton_density, 0.1),  # Ensure positive density
        'V_P': np.maximum(proton_speed, 100),    # Ensure positive speed
        'T_P': np.maximum(proton_temp, 10000)    # Ensure positive temperature
    }

    print("Using simulated Solar Orbiter in-situ and remote sensing data")

###############################################################################
# Create a comprehensive multi-panel plot showing the coordinated observations

fig = plt.figure(figsize=(16, 12))

# Panel 1: EUI image
ax1 = fig.add_subplot(3, 2, 1, projection=eui_map)
norm = ImageNormalize(vmin=eui_map.data.min(), vmax=eui_map.data.max(),
                      stretch=AsinhStretch(0.01))
eui_map.plot(axes=ax1, norm=norm, cmap='sdoaia171')
eui_map.draw_limb(axes=ax1, color='white', linewidth=2)
ax1.set_title(f'Solar Orbiter EUI {eui_map.wavelength}\\n{eui_map.date.strftime("%Y-%m-%d %H:%M")}')

# Panel 2: Magnetic field components
ax2 = fig.add_subplot(3, 2, 2)
ax2.plot(time_obs.datetime, mag_data['B_R'], 'r-', label='B_R', linewidth=1.5)
ax2.plot(time_obs.datetime, mag_data['B_T'], 'g-', label='B_T', linewidth=1.5)
ax2.plot(time_obs.datetime, mag_data['B_N'], 'b-', label='B_N', linewidth=1.5)
ax2.plot(time_obs.datetime, mag_data['B_TOTAL'], 'k-', label='|B|', linewidth=2)
ax2.set_ylabel('Magnetic Field [nT]')
ax2.set_title('MAG: Magnetic Field Components (RTN)')
ax2.legend(loc='upper right')
ax2.grid(True, alpha=0.3)
ax2.xaxis.set_major_formatter(DateFormatter('%H:%M'))

# Panel 3: Proton density
ax3 = fig.add_subplot(3, 2, 3)
ax3.plot(time_obs.datetime, plasma_data['N_P'], 'purple', linewidth=2)
ax3.set_ylabel('Proton Density [cm⁻³]')
ax3.set_title('SWA: Proton Density')
ax3.grid(True, alpha=0.3)
ax3.xaxis.set_major_formatter(DateFormatter('%H:%M'))

# Panel 4: Proton speed
ax4 = fig.add_subplot(3, 2, 4)
ax4.plot(time_obs.datetime, plasma_data['V_P'], 'orange', linewidth=2)
ax4.set_ylabel('Proton Speed [km/s]')
ax4.set_title('SWA: Proton Bulk Speed')
ax4.grid(True, alpha=0.3)
ax4.xaxis.set_major_formatter(DateFormatter('%H:%M'))

# Panel 5: Proton temperature
ax5 = fig.add_subplot(3, 2, 5)
ax5.plot(time_obs.datetime, plasma_data['T_P']/1000, 'brown', linewidth=2)
ax5.set_ylabel('Proton Temperature [1000 K]')
ax5.set_xlabel('Time')
ax5.set_title('SWA: Proton Temperature')
ax5.grid(True, alpha=0.3)
ax5.xaxis.set_major_formatter(DateFormatter('%H:%M'))

# Panel 6: Correlation analysis
ax6 = fig.add_subplot(3, 2, 6)

# Calculate some derived quantities for correlation analysis
dynamic_pressure = 1.67e-6 * plasma_data['N_P'] * (plasma_data['V_P'])**2  # nPa
magnetic_pressure = mag_data['B_TOTAL']**2 / (2 * 4e-7 * np.pi) * 1e-9  # nPa

ax6.scatter(dynamic_pressure, magnetic_pressure, alpha=0.6, s=20)
ax6.set_xlabel('Dynamic Pressure [nPa]')
ax6.set_ylabel('Magnetic Pressure [nPa]')
ax6.set_title('Pressure Balance Analysis')
ax6.grid(True, alpha=0.3)

# Add diagonal line for reference
min_p = min(np.min(dynamic_pressure), np.min(magnetic_pressure))
max_p = max(np.max(dynamic_pressure), np.max(magnetic_pressure))
ax6.plot([min_p, max_p], [min_p, max_p], 'k--', alpha=0.5, label='1:1 line')
ax6.legend()

plt.tight_layout()
plt.show()

###############################################################################
# Perform some basic correlation analysis between in-situ and remote sensing

print("\\n" + "="*60)
print("Coordinated Observation Analysis")
print("="*60)

# Calculate basic statistics
b_mean = np.mean(mag_data['B_TOTAL'])
b_std = np.std(mag_data['B_TOTAL'])
n_mean = np.mean(plasma_data['N_P'])
v_mean = np.mean(plasma_data['V_P'])
t_mean = np.mean(plasma_data['T_P'])

print(f"Magnetic Field Statistics:")
print(f"  Mean |B|: {b_mean:.2f} ± {b_std:.2f} nT")
print(f"  Range: {np.min(mag_data['B_TOTAL']):.2f} - {np.max(mag_data['B_TOTAL']):.2f} nT")

print(f"\\nPlasma Statistics:")
print(f"  Mean proton density: {n_mean:.2f} cm⁻³")
print(f"  Mean proton speed: {v_mean:.1f} km/s")
print(f"  Mean proton temperature: {t_mean/1000:.0f} × 10³ K")

# Look for correlations
correlation_bn = np.corrcoef(mag_data['B_TOTAL'], plasma_data['N_P'])[0,1]
correlation_bv = np.corrcoef(mag_data['B_TOTAL'], plasma_data['V_P'])[0,1]

print(f"\\nCorrelations:")
print(f"  |B| vs Density: {correlation_bn:.3f}")
print(f"  |B| vs Speed: {correlation_bv:.3f}")

###############################################################################
# Highlight unique Solar Orbiter science opportunities

print(f"\\nUnique Solar Orbiter Science Opportunities:")
print(f"{'='*50}")

observer_lon = getattr(eui_map, 'hgln_obs', 45)
observer_lat = getattr(eui_map, 'hglt_obs', 5)
observer_dist = getattr(eui_map, 'dsun_obs', 0.8*u.AU.to(u.m)) / u.AU.to(u.m)

print(f"Observer Position:")
print(f"  Heliographic Longitude: {observer_lon:.1f}°")
print(f"  Heliographic Latitude: {observer_lat:.1f}°")
print(f"  Distance from Sun: {observer_dist:.2f} AU")

print(f"\\nScience Capabilities:")
print(f"  • Simultaneous in-situ and remote sensing observations")
print(f"  • Off-Earth-line perspective for stereoscopic analysis")
print(f"  • Variable heliocentric distance for multi-scale studies")
print(f"  • Direct measurement of solar wind from observed source regions")

# Identify potential science targets
if observer_lon > 30:
    print(f"  • Currently positioned for far-side solar observations")
if abs(observer_lat) > 3:
    print(f"  • Good viewing angle for polar region studies")
if observer_dist < 0.5:
    print(f"  • Close approach enables high-resolution observations")

###############################################################################
# Create a summary plot showing the observational context

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left plot: Orbital context
ax1.set_aspect('equal')
circle = plt.Circle((0, 0), 1, fill=False, linestyle='--', color='blue', alpha=0.5, label='Earth orbit')
ax1.add_patch(circle)

# Plot Sun
ax1.plot(0, 0, 'yo', markersize=15, label='Sun')

# Plot Earth (reference)
ax1.plot(1, 0, 'bo', markersize=8, label='Earth')

# Plot Solar Orbiter
solo_x = observer_dist * np.cos(np.radians(observer_lon))
solo_y = observer_dist * np.sin(np.radians(observer_lon))
ax1.plot(solo_x, solo_y, 'ro', markersize=10, label='Solar Orbiter')

# Draw line of sight
los_length = 2
los_x = [solo_x, solo_x - los_length * np.cos(np.radians(observer_lon))]
los_y = [solo_y, solo_y - los_length * np.sin(np.radians(observer_lon))]
ax1.plot(los_x, los_y, 'r--', alpha=0.7, label='Line of sight')

ax1.set_xlim(-2, 2)
ax1.set_ylim(-2, 2)
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.set_xlabel('Distance [AU]')
ax1.set_ylabel('Distance [AU]')
ax1.set_title('Orbital Configuration')

# Right plot: Data timeline
ax2.plot(time_obs.datetime, mag_data['B_TOTAL'], 'k-', linewidth=2, label='|B| field')
ax2_twin = ax2.twinx()
ax2_twin.plot(time_obs.datetime, plasma_data['N_P'], 'purple', linewidth=2, label='Density')

ax2.set_xlabel('Time')
ax2.set_ylabel('Magnetic Field [nT]', color='black')
ax2_twin.set_ylabel('Density [cm⁻³]', color='purple')
ax2.set_title('In-Situ Measurements Timeline')
ax2.grid(True, alpha=0.3)
ax2.xaxis.set_major_formatter(DateFormatter('%H:%M'))

# Add vertical line for EUI observation time
eui_time = eui_map.date.datetime
ax2.axvline(eui_time, color='red', linestyle=':', linewidth=2, label='EUI image')
ax2.legend(loc='upper left')

plt.tight_layout()
plt.show()

print("\\nThis example demonstrates the power of Solar Orbiter's coordinated")
print("remote sensing and in-situ capabilities for comprehensive solar wind studies.")
