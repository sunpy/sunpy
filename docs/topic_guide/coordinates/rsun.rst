.. _sunpy-topic-guide-coordinates-rsun-in-coordinate-transformations:

**********************************************
The role of rsun in coordinate transformations
**********************************************

In the case of `~sunpy.coordinates.frames.HeliographicCarrington`, one can specify ``observer='self'`` to indicate that the coordinate itself should be used as the observer for defining the coordinate frame.

It is possible to convert from a `~sunpy.coordinates.frames.Helioprojective` frame with one observer location to another `~sunpy.coordinates.frames.Helioprojective` frame with a different observer location.
The transformation requires the coordinate to be 3D, so if the initial coordinate is only 2D, the default assumption maps that 2D coordinate to the surface of the Sun, as defined by the ``rsun`` frame attribute.
The conversion can be performed as follows:

.. code-block:: python

    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord

    >>> from sunpy.coordinates import frames
    >>> import sunpy.coordinates

    >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, observer="earth", obstime="2017-07-26", frame=frames.Helioprojective)
    >>> hpc_out = sunpy.coordinates.Helioprojective(observer="venus", obstime="2017-07-26")
    >>> hpc2 = hpc1.transform_to(hpc_out)

For example if you have two maps, named ``aia`` and ``stereo``:

.. code-block:: python

    >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=aia.coordinate_frame)  # doctest: +SKIP
    >>> hpc2 = hpc1.transform_to(stereo.coordinate_frame)  # doctest: +SKIP
