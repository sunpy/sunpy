"""
=================================
Masking a composite plot
=================================

How to mask out off-disk emission from a composite plot.
"""
# sphinx_gallery_tags = ["Composite", "Masking", "AIA", "Coordinates"]

import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.data.sample
from sunpy.coordinates.utils import coordinate_is_on_solar_disk
from sunpy.map import Map
from sunpy.map.maputils import all_coordinates_from_map

###############################################################################
# Let's import sample data representing the three types of data we want to
# combine into a single plot.

aia = Map(sunpy.data.sample.AIA_171_IMAGE)
hmi = Map(sunpy.data.sample.HMI_MAGNETOGRAM_IMAGE)
euvi = Map(sunpy.data.sample.EUVI_171_IMAGE)

###############################################################################
# We need to reproject the maps to a common grid.  Let's use the AIA map
# as the reference.

common_grid = aia
hmi = hmi.reproject_to(common_grid)
euvi = euvi.reproject_to(common_grid)

###############################################################################
# Now we can create a composite plot.  First, let's create a mask to hide
# off-disk emission.

hpc_coords = all_coordinates_from_map(aia)
mask = coordinate_is_on_solar_disk(hpc_coords)

###############################################################################
# Let's create the composite plot.

fig = plt.figure()
ax = fig.add_subplot(projection=aia)
composite_map = Map(aia.data, aia.meta, mask=mask)
composite_map.plot(axes=ax, title="AIA + HMI + EUVI Composite")
ax.set_xlabel("Helioprojective Longitude [arcsec]")
ax.set_ylabel("Helioprojective Latitude [arcsec]")
plt.colorbar(ax.collections[0], ax=ax, label="Intensity")
plt.show()
