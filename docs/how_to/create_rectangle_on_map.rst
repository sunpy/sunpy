.. _sunpy-how-to-create-rectangle-on-map:

************************************
How to Draw a Rectangle on sunpy map
************************************

``sunpy`` provides a convenient method called :meth:`~sunpy.map.GenericMap.draw_quadrangle` to draw rectangles on maps.
In this guide, we will demonstrate four different methods to draw a rectangle on a `sunpy.map.Map`.

Specify corners with a single `~astropy.coordinates.SkyCoord`
=============================================================

We will use one `~astropy.coordinates.SkyCoord` to represent the two opposite corners.

.. code-block:: python

    >>> import matplotlib.pyplot as plt

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord

    >>> import sunpy.data.sample
    >>> import sunpy.map

    >>> aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE) # doctest: +REMOTE_DATA
    >>> fig = plt.figure(figsize=(5, 5))
    >>> ax = fig.add_subplot(projection=aia_map) # doctest: +SKIP
    >>> aia_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent) # doctest: +SKIP
    >>> coords = SkyCoord(Tx=(100, 500) * u.arcsec, Ty=(200, 500) * u.arcsec,frame=aia_map.coordinate_frame) # doctest: +SKIP
    >>> aia_map.draw_quadrangle(coords, axes=ax, edgecolor="blue") # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

Specify corners with separate `~astropy.coordinates.SkyCoord`
=============================================================

We will use two `~astropy.coordinates.SkyCoord` to represent the two opposite corners.

.. code-block:: python

    >>> aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE) # doctest: +REMOTE_DATA

    >>> fig = plt.figure(figsize=(5, 5))
    >>> ax = fig.add_subplot(projection=aia_map) # doctest: +SKIP
    >>> aia_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)  # doctest: +SKIP

    >>> bottom_left = SkyCoord(-500 * u.arcsec, 200 * u.arcsec, frame=aia_map.coordinate_frame) # doctest: +SKIP
    >>> top_right = SkyCoord(-100 * u.arcsec, 500 * u.arcsec, frame=aia_map.coordinate_frame) # doctest: +SKIP
    >>> aia_map.draw_quadrangle(bottom_left, axes=ax, top_right=top_right, edgecolor="green")  # doctest: +SKIP

    >>> plt.show()  # doctest: +SKIP

Specify one corner with a width and height
==========================================

We will use one `~astropy.coordinates.SkyCoord` to represent the bottom left and supply a width and height to complete the rectangle.

.. code-block:: python

    >>> aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE) # doctest: +REMOTE_DATA

    >>> fig = plt.figure(figsize=(5, 5))
    >>> ax = fig.add_subplot(projection=aia_map) # doctest: +SKIP
    >>> aia_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)  # doctest: +SKIP

    >>> bottom_left = SkyCoord(-500 * u.arcsec, -500 * u.arcsec, frame=aia_map.coordinate_frame) # doctest: +SKIP
    >>> width = 400 * u.arcsec
    >>> height = 300 * u.arcsec
    >>> aia_map.draw_quadrangle(bottom_left, axes=ax, width=width, height=height, edgecolor="yellow")  # doctest: +SKIP

    >>> plt.show()  # doctest: +SKIP

Using pixel coordinates
=======================

We will use a `~astropy.coordinates.SkyCoord` to work out the pixel coordinates instead of using coordinates as we do above.

.. code-block:: python

    >>> aia_map = sunpy.map.Map(sunpy.data.sample.AIA_171_IMAGE) # doctest: +REMOTE_DATA

    >>> fig = plt.figure(figsize=(5, 5))
    >>> ax = fig.add_subplot(projection=aia_map) # doctest: +SKIP
    >>> aia_map.plot(axes=ax, clip_interval=(1, 99.99)*u.percent)  # doctest: +SKIP

    >>> bottom_left = aia_map.wcs.pixel_to_world(600 * u.pixel, 350 * u.pixel) # doctest: +SKIP
    >>> top_right = aia_map.wcs.pixel_to_world(800 * u.pixel, 450 * u.pixel) # doctest: +SKIP
    >>> aia_map.draw_quadrangle(bottom_left, axes=ax, top_right=top_right, edgecolor="red")  # doctest: +SKIP

    >>> plt.show()  # doctest: +SKIP
