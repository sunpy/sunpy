"""
========================
Reading ADAPT FITS Files
========================

This example demonstrates how to load data from the
Air Force Data Assimilative Photospheric Flux Transport (ADAPT) model into a list of `sunpy.map.Map` objects.
"""
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.io import fits

import sunpy.map
from sunpy.coordinates.sun import carrington_rotation_time
from sunpy.net import Fido
from sunpy.net import attrs as a

###############################################################################
# First we will download an ADAPT FITS file.
# To do this we will use the `sunpy.net.Fido` search interface to search for
# ADAPT data for Carrington Rotation 2193, and the first longitude type (0)
# which means the data will be in a Carrington coordinate frame.

date_start = carrington_rotation_time(2193)
date_end = date_start + 6*u.h

result = Fido.search(a.Time(date_start, date_end), a.Instrument('adapt'), a.adapt.ADAPTLonType("0"))
print(result)

downloaded_file = Fido.fetch(result)

###############################################################################
# ADAPT FITS files contain 12 realizations of synoptic magnetogram
# output as a result of varying model assumptions. `This is explained in detail in this talk
# <https://www.swpc.noaa.gov/sites/default/files/images/u33/SWW_2012_Talk_04_27_2012_Arge.pdf>`__.
#
# Because the array in the FITS file is 3D, it cannot be passed directly to `sunpy.map.Map`,
# as this will take the first slice only, ignoring the other realizations.
# Instead, we'll open the FITS file with `astropy.io.fits` and manually pull out each
# model realization.

adapt_fits = fits.open(downloaded_file[0])

###############################################################################
# ``adapt_fits`` is a list of ``HDU`` objects. The first of these contains
# the 12 realizations of the data and a header with sufficient information to build
# a list of maps. We unpack this information into a list of
# ``(data, header)`` tuples where ``data`` are the different adapt realizations.

data_header_pairs = [(map_slice, adapt_fits[0].header) for map_slice in adapt_fits[0].data]

###############################################################################
# Next, pass this list of tuples to `sunpy.map.Map` to make a list of maps:

adapt_maps = sunpy.map.Map(data_header_pairs)

###############################################################################
# ``adapt_maps`` is now a list of our individual adapt realizations.
# Here, we generate a static plot accessing a subset of the individual maps in turn.

fig = plt.figure(figsize=(5, 10), layout='constrained')
for i, a_map in enumerate(adapt_maps[::4]):
    ax = fig.add_subplot(3, 1, i+1, projection=a_map)
    ax.tick_params(axis='x', labelsize=6)
    im = a_map.plot(axes=ax, cmap='bwr', vmin=-50, vmax=50, title=f"Realization {1+i*4:02d}")

fig.colorbar(im, label='Magnetic Field Strength [G]', orientation='horizontal')

plt.show()
