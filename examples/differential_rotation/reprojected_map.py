"""
=============================
Differentially rotating a map
=============================

How to apply differential rotation to a Map.

The example uses the `~sunpy.coordinates.metaframes.RotatedSunFrame` coordinate
metaframe in `sunpy.coordinates` to apply differential rotation to a
Map.  See :ref:`sunpy-coordinates-rotatedsunframe` for more details on
using `~sunpy.coordinates.metaframes.RotatedSunFrame`.

.. note::
   This example requires `reproject` 0.6 or later to be installed.

"""
import matplotlib.pyplot as plt
from reproject import reproject_interp

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

import sunpy.map
from sunpy.coordinates import Helioprojective, RotatedSunFrame, transform_with_sun_center
from sunpy.data.sample import AIA_171_IMAGE

##############################################################################
# First, load an AIA observation.

aiamap = sunpy.map.Map(AIA_171_IMAGE)
in_time = aiamap.date

##############################################################################
# Let's define the output frame to be five days in the future for an observer
# at Earth (i.e., about five degrees offset in heliographic longitude compared
# to the location of AIA in the original observation).

out_time = in_time + 5*u.day
out_frame = Helioprojective(observer='earth', obstime=out_time)

##############################################################################
# For the reprojection, the definition of the target frame can be
# counter-intuitive. Reprojection should be performed between two frames that
# point to inertial locations at the same time, but ``out_frame`` is not at
# the same time as the original frame (``out_time`` versus ``in_time``).
# The `~sunpy.coordinates.metaframes.RotatedSunFrame` metaframe allows one to
# specify coordinates in a coordinate frame at one time (here, ``out_frame``)
# to refer to inertial locations at a different time (here, ``in_time``,
# which is supplied through the keyword argument ``rotated_time``).
# `~sunpy.coordinates.metaframes.RotatedSunFrame` will account for the
# appropriate amount of solar (differential) rotation for the time
# difference between ``out_time`` and ``in_time``.

rot_frame = RotatedSunFrame(base=out_frame, rotated_time=in_time)
print(rot_frame)

##############################################################################
# Construct a WCS object for the output map with the target
# ``RotatedSunHelioprojective`` frame specified instead of the regular
# `~sunpy.coordinates.frames.Helioprojective` frame.
# If one has an actual ``Map`` object at the desired output
# time (e.g., the actual AIA observation at the output time), one can use the
# WCS object from that ``Map`` object (e.g., ``mymap.wcs``) instead of
# constructing a custom WCS.

out_shape = aiamap.data.shape
out_center = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=out_frame)
header = sunpy.map.make_fitswcs_header(out_shape,
                                       out_center,
                                       scale=u.Quantity(aiamap.scale))
out_wcs = WCS(header)
out_wcs.coordinate_frame = rot_frame

##############################################################################
# Reproject the map from the input frame to the output frame.  We use the
# :func:`~sunpy.coordinates.transform_with_sun_center` context manager so that
# the coordinates stay defined relative to Sun center as the Sun moves slowly
# over time.

with transform_with_sun_center():
    arr, _ = reproject_interp(aiamap, out_wcs, out_shape)

##############################################################################
# Finally, we create the output map and preserve the original map's plot
# settings.

out_warp = sunpy.map.Map(arr, out_wcs)
out_warp.plot_settings = aiamap.plot_settings

##############################################################################
# Let's plot the differentially rotated Map next to the original Map.

fig = plt.figure(figsize=(12, 4))

ax1 = fig.add_subplot(1, 2, 1, projection=aiamap)
aiamap.plot(vmin=0, vmax=20000, title='Original map')
plt.colorbar()

ax2 = fig.add_subplot(1, 2, 2, projection=out_warp)
out_warp.plot(vmin=0, vmax=20000,
              title=f"Reprojected to an Earth observer {(out_time - in_time).to('day')} later")
plt.colorbar()

plt.show()
