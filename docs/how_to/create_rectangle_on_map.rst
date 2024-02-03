.. _sunpy-how-to-create-rectangle-on-map:

************************************
How to Draw a Rectangle on sunpy map
************************************


SunPy provides a convenient method called :meth:`~sunpy.map.GenericMap.draw_quadrangle` to draw rectangles on maps. 
In this guide, we'll demonstrate four different ways to draw a rectangle on a map using SunPy.

Specify Corners as a Single `~astropy.coordinates.SkyCoord` Object
==================================================================

Taking two opposite corners as a single SkyCoord object

.. code-block:: python

    >>> import matplotlib.pyplot as plt

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord

    >>> import sunpy.data.sample
    >>> import sunpy.map

    >>> aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> fig = plt.figure(figsize=(5, 5))
    >>> ax = fig.add_subplot(projection=aia_map)
    >>> aia_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
    >>> coords = SkyCoord(
        Tx=(100, 500) * u.arcsec,
        Ty=(200, 500) * u.arcsec,
        frame=aia_map.coordinate_frame,
        )
    >>> aia_map.draw_quadrangle(
        coords,
        axes=ax,
        edgecolor="blue",
        linestyle="-",
        linewidth=2,
        label='2-element SkyCoord'
        )

Specify Corners as Separate `~astropy.coordinates.SkyCoord` Objects
===================================================================

Taking two opposite corners as separate SkyCoord Objects

.. code-block:: python

    >>> aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> fig = plt.figure(figsize=(5, 5))
    >>> ax = fig.add_subplot(projection=aia_map)
    >>> aia_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
    >>> bottom_left = SkyCoord(-500 * u.arcsec, 200 * u.arcsec, frame=aia_map.coordinate_frame)
    >>> top_right = SkyCoord(-100 * u.arcsec, 500 * u.arcsec, frame=aia_map.coordinate_frame)
    >>> aia_map.draw_quadrangle(
        bottom_left,
        axes=ax,
        top_right=top_right,
        edgecolor="green",
        linestyle="--",
        linewidth=2,
        label='two SkyCoords'
        )

Specify One Corner, Width, and Height
=====================================

.. code-block:: python

    >>> aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> fig = plt.figure(figsize=(5, 5))
    >>> ax = fig.add_subplot(projection=aia_map)
    >>> aia_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
    >>> bottom_left = SkyCoord(-500 * u.arcsec, -500 * u.arcsec, frame=aia_map.coordinate_frame)
    >>> width = 400 * u.arcsec
    >>> height = 300 * u.arcsec
    >>> aia_map.draw_quadrangle(
        bottom_left,
        axes=ax,
        width=width,
        height=height,
        edgecolor="yellow",
        linestyle="-.",
        linewidth=2,
        label='width/height'
        )

Draw a Rectangle in Pixel Coordinates
=====================================

.. code-block:: python

    >>> aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE)
    >>> fig = plt.figure(figsize=(5, 5))
    >>> ax = fig.add_subplot(projection=aia_map)
    >>> aia_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)
    >>> bottom_left = aia_map.wcs.pixel_to_world(600 * u.pixel, 350 * u.pixel)
    >>> top_right = aia_map.wcs.pixel_to_world(800 * u.pixel, 450 * u.pixel)
    >>> aia_map.draw_quadrangle(
        bottom_left,
        axes=ax,
        top_right=top_right,
        edgecolor="red",
        linestyle=":",
        linewidth=2,
        label='pixel_to_world()'
        )