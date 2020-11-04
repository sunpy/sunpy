"""
========================================
Saving and loading coordinates with asdf
========================================

In this example we are going to look at saving and loading collections of
coordinates with `asdf <https://asdf.readthedocs.io/en/latest/>`__.

asdf is a modern file format designed to meet the needs of the astronomy
community. It has deep integration with Python and SunPy and Astropy as well as
implementations in other languages. It can be used to store known Python
objects in a portable, well defined file format. It is primarily useful for
storing complex Astropy and SunPy objects in a way that can be loaded back into
the same form as they were saved.

.. note::
    This example requires Astropy 3.2 and asdf 2.3.0

"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

import asdf
import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.coordinates import frames
from sunpy.data.sample import AIA_171_IMAGE
from sunpy.sun import constants

################################################################################
# To get started let's use a function to get the coordinates of a semi-circular
# loop from
# `this <https://sunpy.org/posts/2018/2018-07-21-coronal-loop-coordinates.html>`__
# blog post by Will Barnes to generate ourselves some coordinates.


@u.quantity_input
def semi_circular_loop(length: u.m, latitude: u.deg = 0*u.deg):
    """
    Return a Heliographic Stonyhurst coordinate object with points of a semi circular loop in it.
    """
    r_sun = constants.radius

    def r_2_func(x):
        return np.arccos(0.5 * x / r_sun.to(u.cm).value) - np.pi + length.to(u.cm).value / 2. / x

    # Find the loop radius corresponding to the loop length
    r_2 = scipy.optimize.bisect(r_2_func,
                                length.to(u.cm).value / (2 * np.pi),
                                length.to(u.cm).value / np.pi) * u.cm
    alpha = np.arccos(0.5 * (r_2 / r_sun))
    phi = np.linspace(-np.pi * u.rad + alpha, np.pi * u.rad - alpha, 2000)

    hcc_frame = frames.Heliocentric(
        observer=frames.HeliographicStonyhurst(lon=0 * u.deg, lat=latitude, radius=1 * u.AU))

    return SkyCoord(
        x=r_2 * np.sin(phi),
        y=0 * u.cm,
        z=r_2 * np.cos(phi) + r_sun,
        frame=hcc_frame).transform_to('heliographic_stonyhurst')


################################################################################
# Use this function to generate a `~astropy.coordinates.SkyCoord` object.
loop_coords = semi_circular_loop(500*u.Mm, 30*u.deg)
print(loop_coords.shape)
# print the first and last coordinate point
print(loop_coords[[0, -1]])


################################################################################
# This is a regular coordinate object that can be transformed to other frames
# or overplotted on images. For instance we could overplot it on an AIA image

aiamap = sunpy.map.Map(AIA_171_IMAGE)

ax = plt.subplot(projection=aiamap)
aiamap.plot(axes=ax, clip_interval=(1, 99.5) * u.percent)
ax.plot_coord(loop_coords, 'r')

plt.show()


################################################################################
# We can now save these loop points to an asdf file to use later. The advantage
# of saving them to asdf is that all the metadata about the coordinates will be
# preserved, and when we load the asdf, we will get back an identical
# `~astropy.coordinates.SkyCoord` object.
#
# asdf files save a dictionary to a file, so to save the loop coordinates we
# need to put them into a dictionary. This becomes what asdf calls a tree.

tree = {'loop_points': loop_coords}

with asdf.AsdfFile(tree) as asdf_file:
    asdf_file.write_to("loop_coords.asdf")


################################################################################
# This asdf file is a portable file and can be safely loaded by anyone with
# Astropy and SunPy installed. We can reload the file like so:

with asdf.open("loop_coords.asdf") as input_asdf:
    new_coords = input_asdf['loop_points']

print(new_coords.shape)
# print the first and last coordinate point
print(new_coords[[0, -1]])
