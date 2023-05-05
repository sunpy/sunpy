.. _topic-guide-rsun:

*********************************
The role of ``rsun`` in ``sunpy``
*********************************

When constructing a ``map.wcs`` a key component is ``rsun``, the height at which the emission in the data is presumed to originate.
Data providers often provide this information to end users via the fits keyword: ``RSUN_REF``,
where this is not available ``sunpy`` assumes a standard value defined by `sunpy.sun.constants.radius`.

The intention behind ``rsun`` is to define the height at which the emission originates and therefore provide 
3-dimensional information about the location of the detection. In practice, the radiation detected by a single filter,
especially for coronal emission, comes from a variety of heights and this range of heights is not necessarily represented by ``rsun``.

Reprojecting between frames with different ``rsun``
===================================================

When the values of ``rsun`` for two wcs instances 


Transforming between ``helioprojective`` frames
===============================================

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

An example with two maps, named ``aia`` and ``stereo``::

.. code-block:: python
  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=aia.coordinate_frame)  # doctest: +SKIP
  >>> hpc2 = hpc1.transform_to(stereo.coordinate_frame)  # doctest: +SKIP