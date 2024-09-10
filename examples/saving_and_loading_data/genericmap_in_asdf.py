"""
=======================================
Saving and loading sunpy Maps with ASDF
=======================================

In this example, we are going to look at how we can save and load a
`~sunpy.map.GenericMap` with `asdf <https://asdf.readthedocs.io/en/latest/>`__.

ASDF is a modern file format designed to meet the needs of the astronomy
community. It has deep integration with Python, SunPy, and Astropy, as well as
implementations in other languages. It can be used to store known Python
objects in a portable, well-defined file format. It is primarily useful for
storing complex Astropy and SunPy objects in a way that can be loaded back into
the same form as they were saved.

Here, even though we will be working with `~sunpy.map.sources.sdo.AIAMap`
specifically, the process can be extended to any `~sunpy.map.GenericMap`,
including ones created using custom FITS files.
"""
import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.data.sample
import sunpy.map

################################################################################
# We begin by creating an `~sunpy.map.sources.sdo.AIAMap` object using the
# sample data.

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
aia_map.peek(clip_interval=(1, 99.99) * u.percent)

################################################################################
# We can now save this object to an ASDF file to use later. Saving it like this
# allows us to preserve all of the metadata of the object along with the actual
# array data. When we load the ASDF file again, we get an identical
# `~sunpy.map.sources.sdo.AIAMap` object.

# Save the map to an ASDF file
aia_map.save("sunpy_map.asdf")

################################################################################
# This ASDF file is a portable file and can be safely loaded by anyone with
# astropy, sunpy, and asdf installed. We can reload it using the `sunpy.map`:

# Load the ASDF file
pairs = sunpy.map.Map("sunpy_map.asdf")

pairs.plot(clip_interval=(1, 99.99) * u.percent)

plt.show()
