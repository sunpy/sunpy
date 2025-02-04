"""
===========================================
Time-Distance Slice Plot of Sequence of Map
===========================================

This example showcases how you can sunpy.map.pixelate_coord_path() and sunpy.map.sample_at_coords() on a sequence of images to extract a time-distance slice to identify features in solar atmosphere
"""
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, SqrtStretch

import sunpy.map
from sunpy.coordinates import propagate_with_solar_surface
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# First we will import sequence of images from AIA 171 angstrom separated by 1 hour,
# and create a sequence of maps.
query = Fido.search(a.Time('2018-05-30 02:00:00', '2018-05-30 10:00:00'),
                    a.Instrument.aia,
                    a.Wavelength(171*u.angstrom),
                    a.Sample(1*u.hr))
print(query)
files = Fido.fetch(query)
aia_sequence = sunpy.map.Map(files, sequence=True)

########################################################################################################################
# sunpy.map.pixelate_coord_path will return the pixel coordinates for every pixel that intersects with a coordinate path.
# To obtain the values of these pixels/intensity, pass the output to sunpy.map.sample_at_coords().

line_coords = SkyCoord([-350, 100], [-250, 650], unit=u.arcsec, frame=aia_sequence[0].coordinate_frame)
intensity_coords = sunpy.map.pixelate_coord_path(aia_sequence[0], line_coords)
intensity = sunpy.map.sample_at_coords(aia_sequence[0], intensity_coords)
angular_separation = intensity_coords.separation(intensity_coords[0]).to(u.arcsec)


###############################################################################
# Let's now create rectangular slice i.e map cutout.

corner = SkyCoord(Tx=-250*u.arcsec, Ty=0*u.arcsec, frame=aia_sequence[0].coordinate_frame)
width = 250*u.arcsec
cutout_map = aia_sequence[0].submap(corner, width=width, height=500*u.arcsec)

##################################################################################################
# Now that we have our cutout, we can reproject each map in our sequence to
# the WCS of that cutout, using the `~sunpy.coordinates.propagate_with_solar_surface`
# context manager to adjust the field of view of the cutout with the rotation of the solar surface.

with propagate_with_solar_surface():
    aia_sequence_aligned = sunpy.map.Map([m.reproject_to(cutout_map.wcs) for m in aia_sequence], sequence=True)

################################################################################################
#Now we can plot the time-distance slice of the sequence of maps using matplotlib.pyplot.imshow().


fig, ax = plt.subplots(figsize=(15, 6))
norm = norm=ImageNormalize(vmin=0, vmax=3e3, stretch=SqrtStretch())

x_offset = 0
time_ticks = []
time_labels = []

for i, m in enumerate(aia_sequence_aligned):
    extent = [x_offset, x_offset + m.data.shape[1], 0, m.data.shape[0]]
    ax.imshow(m.data, origin='lower', cmap='sdoaia171', extent=extent, norm=norm)

    ax.plot(intensity_coords.Tx.value - corner.Tx.value + x_offset, intensity_coords.Ty.value - corner.Ty.value, color='white')

    time_ticks.append(x_offset + m.data.shape[1] / 2)
    time_labels.append(aia_sequence[i].date.strftime('%H:%M'))

    x_offset += m.data.shape[1]

ax.set_xlim(0, x_offset)
ax.set_ylim(0, m.data.shape[0])

ax.set_xlabel('Time')
ax.set_ylabel('Arcsec along slit')
ax.set_title('Time distance along slit')

ax.set_xticks(time_ticks)
ax.set_xticklabels(time_labels)

plt.tight_layout()
plt.show()
