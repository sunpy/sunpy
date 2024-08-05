"""
Parsing ADAPT Ensemble FITS file
================================

This example demonstrates how to parse an ADAPT FITS file into a `sunpy.map.MapSequence`.
"""
import matplotlib.pyplot as plt
from matplotlib import gridspec

import astropy.units as u
from astropy.io import fits

import sunpy.map
from sunpy.coordinates.sun import carrington_rotation_time
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# First we will acquire a ADAPT fits file.
#
# To do this we will use the `sunpy.net.Fido` search interface to search for
# ADAPT data. We will search for ADAPT data for Carrington Rotation 2193, and
# request the ADAPT data for the first longitude type (0). This will return a
# list of results, each of which is a realization of the ADAPT data. We will
# then download the first result.

date_start = carrington_rotation_time(2193)
date_end = date_start + (3*2 * u.hour)
print(date_start, date_end)
result = Fido.search(a.Time(date_start, date_end), a.Instrument('adapt'), a.adapt.ADAPTLonType("0"))
downloaded_file = Fido.fetch(result)

###############################################################################
# ADAPT synoptic magnetograms contain 12 realizations of synoptic magnetograms
# output as a result of varying model assumptions. `This is explained in detail by a talk
# from 2012 <https://www.swpc.noaa.gov/sites/default/files/images/u33/SWW_2012_Talk_04_27_2012_Arge.pdf>`__.
#
# Because the FITS data is 3D, it cannot be passed directly to `sunpy.map.Map`,
# because this will take the first slice only and the other realizations are
# lost. We want to end up with a `sunpy.map.MapSequence` containing all these
# realizations as individual maps. These maps can then be individually accessed
# and PFSS solutions generated from them.
#
# Now we will open the FITS file:

adapt_fits = fits.open(downloaded_file[0])

###############################################################################
# ``adapt_fits`` is a list of ``HDPair`` objects. The first of these contains
# the 12 realizations data and a header with sufficient information to build
# the `~sunpy.map.MapSequence`. We unpack this ``HDPair`` into a list of
# ``(data, header)`` tuples where ``data`` are the different adapt realizations.

data_header_pairs = [(map_slice, adapt_fits[0].header) for map_slice in adapt_fits[0].data]

###############################################################################
# Next, pass this list of tuples as the argument to `sunpy.map.Map` to create
# the map sequence:

adapt_sequence = sunpy.map.Map(data_header_pairs, sequence=True)

###############################################################################
# ``adapt_map_sequence`` is now a list of our individual adapt realizations.
# Here, we generate a static plot accessing the individual maps in turn:

fig = plt.figure(figsize=(20, 15))
gs = gridspec.GridSpec(4, 3, figure=fig)
for i, a_map in enumerate(adapt_sequence):
    ax = fig.add_subplot(gs[i], projection=a_map)
    ax.tick_params(axis='x', labelsize=6)
    a_map.plot(axes=ax, cmap='bwr', vmin=-2, vmax=2, title=f"Realization {1+i:02d}")

fig.tight_layout()

plt.show()
