"""
=========================
Map metadata modification
=========================

How to query map metadata for changes.

sunpy has a series of map sources that can fix common issues with the original metadata
stored in FITS files. A copy of the original (unaltered) metadata is stored, so
any changes that sunpy (or the user) subsequently makes to the metadata can be
easily queried.

In the example below, we load a HMI sample image, and query the metadata for
any added, removed, or modified items.
"""
import astropy.units as u

import sunpy.map
from sunpy.data.sample import HMI_LOS_IMAGE

###############################################################################
# Start by creating a map.

hmimap = sunpy.map.Map(HMI_LOS_IMAGE)

###############################################################################
# Now query the ``.meta`` attribute for any changes. We can see that nothing
# has been added or removed, but the 'bunit' key has been updated from "Gauss"
# to "G" to make it FITS standard compliant.

print("Added items:", hmimap.meta.added_items)
print("Removed items:", hmimap.meta.removed_items)
print("Modified items:", hmimap.meta.modified_items)

###############################################################################
# If we modify the map in a way that updates the metadata, these properties
# allow use to easily see what has changed. As an example, lets rotate the map
# by 90 degrees, and see what has been updated.

hmimap = hmimap.rotate(90 * u.deg)
print("Added items:", hmimap.meta.added_items)
print("Removed items:", hmimap.meta.removed_items)
print("Modified items:", hmimap.meta.modified_items)
