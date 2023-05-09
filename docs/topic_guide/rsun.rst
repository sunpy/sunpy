.. _topic-guide-rsun:

*********************************
The role of ``rsun`` in ``sunpy``
*********************************

When a ``map.wcs`` is constructed, a key variable is ``rsun``; the height at which emission originates.
Data providers often include this information for end users via the fits keyword: ``RSUN_REF``,
where this keyword is not available ``sunpy`` assumes a standard value defined by `sunpy.sun.constants`.

In theory, the radiation detected by a single filter, especially for coronal emission,
comes from a variety of heights and ``rsun`` is only a loose approximation.
In practice though, the value of ``rsun`` is important when working with coordinates and transforms, it provides 3 dimensional information about the location of a map pixel.
Idiosyncrasies can occur which can be pre-empted or fixed by an understanding of behaviors surrounding ``rsun``.

Reprojecting between frames with different ``rsun``
===================================================

When the values of ``rsun`` from two wcs instances are different, issues with reprojecting between those frames can be encountered:
:meth:`~sunpy.map.GenericMap.reproject_to` by default enforces a round-trip behavior,
the idea being that you should only trust the reprojection if the 2D coordinates from each observer both resolve to the same 3D point (within a pixel volume).
When the values for ``rsun`` are different, that criterion fails towards the limb.
In other cases, this furthers results in banding as the criterion fails, then succeeds and then fails again.

Changing the value of ``rsun`` can fix this issue:
Take the example, :ref:`sphx_glr_generated_gallery_map_transformations_reprojection_different_observers.py`
If we run the reprojection twice, once before fixing the discrepancy in the ``rsun_ref`` metadata and once after:

.. code-block:: python

    >>> import sunpy.map
    >>> from sunpy.map import Map
    >>> from sunpy.data.sample import AIA_193_JUN2012, STEREO_A_195_JUN2012
    >>> map_aia = Map(AIA_193_JUN2012)
    >>> map_euvi = Map(STEREO_A_195_JUN2012)
    >>> outmap1 = map_euvi.reproject_to(map_aia.wcs)
    >>> map_euvi.meta['rsun_ref'] = map_aia.meta['rsun_ref']
    >>> outmap2 = map_euvi.reproject_to(map_aia.wcs)

we can see the difference in the appearance of the reprojected maps at the limb.

Transforming between ``helioprojective`` frames
===============================================

It is possible to convert from a `~sunpy.coordinates.frames.Helioprojective` frame with one observer location to another `~sunpy.coordinates.frames.Helioprojective` frame with a different observer location.
The transformation requires the coordinate to be 3D, so if the initial coordinate is only 2D, the default assumption maps that 2D coordinate to the surface of the Sun, at a radius defined by the ``rsun`` frame attribute.
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

  >>> hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=aia.coordinate_frame)  # doctest: +SKIP
  >>> hpc2 = hpc1.transform_to(stereo.coordinate_frame)  # doctest: +SKIP
