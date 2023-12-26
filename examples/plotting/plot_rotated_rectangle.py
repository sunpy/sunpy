import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
import sunpy.data.sample
import sunpy.map

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)

rotation_angle = 90 * u.deg
center_coord = SkyCoord(0 * u.arcsec, 0 * u.arcsec, frame=aia_map.coordinate_frame)

width = 400 * u.arcsec
height = 300 * u.arcsec

offset_frame = SkyOffsetFrame(origin=center_coord, rotation=rotation_angle)

corner_bottom_left = SkyCoord(lon=-width / 2, lat=-height / 2, frame=offset_frame)
corner_top_right = SkyCoord(lon=width / 2, lat=height / 2, frame=offset_frame)

corner_bottom_left = corner_bottom_left.transform_to(aia_map.coordinate_frame)
corner_top_right = corner_top_right.transform_to(aia_map.coordinate_frame)

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(projection=aia_map)
aia_map.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)

aia_map.draw_quadrangle(
    bottom_left=corner_bottom_left,
    top_right=corner_top_right,
    axes=ax,
    edgecolor="purple",
    linestyle="--",
    linewidth=2,
    label='Rotated Rectangle'
)

ax.legend()
plt.show()