"""
=====================
Artemis II trajectory
=====================

This example visualizes the trajectory of the Artemis II spacecraft.

Artemis II was NASA's first crewed mission to the Moon since Apollo 17 in 1972,
flying four astronauts around the Moon on a ~10 day test flight in April 2026.
The trajectory is visualized in two different coordinate frames.  The plots
also highlight the segment of the trajectory when Artemis II was in eclipse.
"""
# sphinx_gallery_tags = ["Coordinates", "Solar Eclipse"]

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter

import astropy.units as u
from astropy.constants import R_earth
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time

from sunpy.coordinates import get_horizons_coord, sun
from sunpy.time import parse_time

##############################################################################
# First, define times spanning the Artemis II mission, with higher resolution
# across eclipse transitions (i.e., the four contacts).  The contact times are
# approximate.

t_start = parse_time("2026-Apr-02 01:58:33")
t_c1 = parse_time("2026-Apr-07 00:32:46")  # start of partial eclipse
t_c2 = parse_time("2026-Apr-07 00:34:28")  # start of total eclipse
t_c3 = parse_time("2026-Apr-07 01:28:55")  # end of total eclipse
t_c4 = parse_time("2026-Apr-07 01:31:03")  # end of partial eclipse
t_end = parse_time("2026-Apr-10 23:54:22")

time_spans = [(t_start, t_c1 - 10*u.s, 5*u.min),
              (t_c1 - 30*u.s, t_c2 + 30*u.s, 5*u.s),
              (t_c2 + 30*u.s, t_c3 - 30*u.s, 5*u.min),
              (t_c3 - 30*u.s, t_c4 + 30*u.s, 5*u.s),
              (t_c4 + 30*u.s, t_end, 5*u.min)]
times = Time(np.concatenate([np.arange(t1.jd, t2.jd, dt.to_value(u.day))
                             for t1, t2, dt in time_spans]), format='jd')

##############################################################################
# Use JPL Horizons via :func:`~sunpy.coordinates.get_horizons_coord` to
# retrieve relevant coordinates, and then use the coordinate framework to
# convert to ecliptic coordinates.  There is no convenient coordinate frame
# for the Earth-Moon system, so convert to a heliocentric frame for now.

# Use a JPL ephemeris because astropy's built-in ephemeris is not accurate enough
solar_system_ephemeris.set('de440s')

# Get the Artemis II coordinate in Heliocentric Mean Ecliptic
artemis = get_horizons_coord("Artemis II", times).heliocentricmeanecliptic

# Get other relevant coordinates
# Specify NAIF IDs for Earth and Moon due to multiple matches for string input
earth = get_horizons_coord(399, times).heliocentricmeanecliptic
moon = get_horizons_coord(301, times).heliocentricmeanecliptic
em_barycenter = get_horizons_coord("Earth-Moon Barycenter", times).heliocentricmeanecliptic

##############################################################################
# Determine when Artemis II was in eclipse using :func:`sunpy.coordinates.sun.eclipse_amount`.
# When the eclipse percentage is greater than 0, at least part of the Sun is
# eclipsed by the Moon.  Be aware that this function assumes a uniform lunar
# radius, but features of the lunar terrain may be comparable to the apparent
# size of the Sun as seen from Artemis II, so the calculation is only an
# approximation.

eclipse_percentage = sun.eclipse_amount(artemis)
eclipsed = np.flatnonzero(eclipse_percentage > 0)  # at least partially eclipsed

##############################################################################
# Plot the eclipse percentage when transitioning in an out of eclipse.

fig, axs = plt.subplots(1, 2, layout="constrained")

axs[0].plot(times.datetime64, eclipse_percentage, '.-')
axs[0].set_xlim((t_c1 - 1*u.min).datetime64, (t_c2 + 1*u.min).datetime64)
axs[0].set_title("Entering eclipse")

axs[1].plot(times.datetime64, eclipse_percentage, '.-')
axs[1].set_xlim((t_c3 - 1*u.min).datetime64, (t_c4 + 1*u.min).datetime64)
axs[1].set_title("Exiting eclipse")

for ax in axs:
    ax.grid()
    ax.xaxis.set_major_formatter(DateFormatter('%m-%d %H:%M'))
    ax.tick_params('x', rotation=90)
    ax.set_ylabel("Eclipse percentage")

##############################################################################
# Shift the coordinates to have the Earth-Moon barycenter as the origin,
# convert to units of Earth radii, and keep only the X and Y components for
# later plotting.

artemis_x, artemis_y, _ = ((artemis.cartesian - em_barycenter.cartesian).xyz / R_earth).decompose()
earth_x, earth_y, _ = ((earth.cartesian - em_barycenter.cartesian).xyz / R_earth).decompose()
moon_x, moon_y, _ = ((moon.cartesian - em_barycenter.cartesian).xyz / R_earth).decompose()

##############################################################################
# Plot the Artemis II trajectory in fixed ecliptic X-Y coordinates. The motion
# of the Earth relative to the Earth-Moon barycenter is not discernible on
# this plot.  The segment of the trajectory when Artemis II was in eclipse is
# highlighted.

fig, ax = plt.subplots()

ax.plot(earth_x, earth_y, ls='dashed', color='b', label='Earth')
ax.plot(earth_x[-1], earth_y[-1], '.', color='b')

ax.plot(moon_x, moon_y, ls='dashed', color='g', label='Moon')
ax.plot(moon_x[-1], moon_y[-1], '.', color='g')

ax.plot(artemis_x, artemis_y, color='k', label='Artemis II')
ax.plot(artemis_x[-1], artemis_y[-1], '.', color='k')
ax.plot(artemis_x[eclipsed], artemis_y[eclipsed], color='m', lw=3, label='eclipse')

ax.set_title('Fixed coordinate frame')
ax.set_xlabel('Ecliptic X (Earth radii)')
ax.set_ylabel('Ecliptic Y (Earth radii)')
ax.set_aspect('equal')
ax.legend(loc='center right')

##############################################################################
# Transform the X and Y components so that the we are in the frame co-rotating
# with the Moon's orbital motion.

angle = np.arctan2(moon_y, moon_x)
c, s = np.cos(-angle), np.sin(-angle)

artemis_xp, artemis_yp = artemis_x * c - artemis_y * s, artemis_x * s + artemis_y * c
earth_xp, earth_yp = earth_x * c - earth_y * s, earth_x * s + earth_y * c
moon_xp, moon_yp = moon_x * c - moon_y * s, moon_x * s + moon_y * c

##############################################################################
# Plot the Artemis II trajectory in coordinates co-rotating with the Moon's
# orbital motion.

fig, ax = plt.subplots()

ax.plot(earth_xp, earth_yp, ls='dashed', color='b', label='Earth')
ax.plot(earth_xp[-1], earth_yp[-1], '.', color='b')

ax.plot(moon_xp, moon_yp, ls='dashed', color='g', label='Moon')
ax.plot(moon_xp[-1], moon_yp[-1], '.', color='g')

ax.plot(artemis_xp, artemis_yp, color='k', label='Artemis II')
ax.plot(artemis_xp[-1], artemis_yp[-1], '.', color='k')
ax.plot(artemis_xp[eclipsed], artemis_yp[eclipsed], color='m', lw=3, label='eclipse')

ax.set_title('Coordinate frame co-rotating with the Moon')
ax.set_xlabel('X (Earth radii)')
ax.set_ylabel('Y (Earth radii)')
ax.set_aspect('equal')
ax.legend(loc='center')

plt.show()

# sphinx_gallery_thumbnail_number = 3
