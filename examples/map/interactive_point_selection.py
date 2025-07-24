"""
=============================================
Interactively Selecting Coordinates on a Plot
=============================================

In this example we demonstrate how to extract world coordinates of points selected interactively.

.. note::

    We can not demonstrate the interactivity in the documentation.
    To see this example work correctly please download it and run it locally.

"""

import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.data.sample
import sunpy.map

###############################################################################
# We will start with using sunpy's sample data for this example.

aia_map = sunpy.map.Map(sunpy.data.sample.AIA_193_IMAGE)

fig = plt.figure()
ax = fig.add_subplot(projection=aia_map)
aia_map.plot(axes=ax)


###############################################################################
# Next we set up a callback function that matplotlib will call each
# time a button is pressed.  We are going to store the clicked
# positions in world coordinates to make them useful for plotting on
# other maps or extracting timeseries.

clicked_points = []


def on_click(event):
    # We only want to react to left mouse clicks inside the image
    if event.button is MouseButton.LEFT and event.inaxes:
        # The xdata,ydata point is the pixel coordinates of the displayed map
        point = (event.xdata, event.ydata)

        # We can plot the point
        plt.plot(*point, "o")
        # We have to redraw to update the plot
        plt.draw()

        # Convert the pixel point to world coordinates and store them
        clicked_points.append(aia_map.wcs.pixel_to_world(*point))


# Connect the event to the plot
plt.connect("button_press_event", on_click)
# Display the figure
plt.show()


###############################################################################
# We are now going to generate an example list of points for the documentation.
# You can remove this code if running locally (or ignore it)
if not clicked_points:
    clicked_points = SkyCoord(
        [-270.75768304, -163.4167166, -78.98673302, 24.88601012, 74.19676398, 108.62529128, 173.31041542],
        [213.16666796, 202.67939786, 205.31559493, 230.42887498, 247.26223829, 298.72639623, 317.47242019],
        unit=u.arcsec,
        frame=aia_map.coordinate_frame,
    )

###############################################################################
# If points were selected on the plot then we have a list of
# `~astropy.coordinates.SkyCoord` objects.
# So we stack the list of clicked_points into a single SkyCoord

clicked_points = SkyCoord(clicked_points)
print(clicked_points)

###############################################################################
# Now pixelate the path and extract the intensities
pixel_path = sunpy.map.pixelate_coord_path(aia_map, clicked_points)
intensities = sunpy.map.sample_at_coords(aia_map, pixel_path)

fig = plt.figure()
ax = fig.add_subplot()
ax.plot(intensities)
ax.set_xlabel("pixel along path")
ax.set_ylabel("Intensity along path")
plt.show()
